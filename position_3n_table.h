/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-3N.
 *
 * HISAT-3N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-3N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef POSITION_3N_TABLE_H
#define POSITION_3N_TABLE_H

#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <thread>
#include <cassert>
#include "alignment_3n_table.h"

using namespace std;

extern bool CG_only;
extern long long int loadingBlockSize;

/**
 * store unique information for one base information with readID, and the quality.
 */
class uniqueID
{
public:
    unsigned long long readNameID;
    bool isConverted;
    char quality;
    bool removed;

    uniqueID(unsigned long long InReadNameID,
             bool InIsConverted,
             char& InQual){
        readNameID = InReadNameID;
        isConverted = InIsConverted;
        quality = InQual;
        removed = false;
    }
};

/**
 * basic class to store reference position information
 */
class Position{
    mutex mutex_;
public:
    short chromosomeId;
    long long int location; // 1-based position
    char strand; // +(REF) or -(REF-RC)
    int convertedCount = 0;
    int unconvertedCount = 0;
    bool empty = true;

    void initialize() {
        chromosomeId = -1;
        location = -1;
        strand = '?';
        convertedCount = 0;
        unconvertedCount = 0;
        empty = true;
    }

    Position(){
        initialize();
    };

    /**
     * return true if there is mapping information in this reference position.
     */
    bool isEmpty() const {
        return empty;
    }

    /**
     * set the chromosome, location (position), and strand information.
     */

    void set (int InChromosomeId, long long int inputLoc) {
        chromosomeId = InChromosomeId;
        location = inputLoc + 1;
    }

    void set(char inputStrand) {
        strand = inputStrand;
    }

    /**
     * append the SAM information into this position.
     */
    void appendBase (PosQuality& input, Alignment& a) {
            if (empty) empty = false;
            if (input.converted) {
                convertedCount++;
            } else {
                unconvertedCount++;
            }
    }

};

/**
 * store all reference position in this class.
 */
class Positions{
public:
    vector<Position*> refPositions; // the pool of all current reference position.
    string chromosome; // current reference chromosome name.'
    int curChromosomeId;
    long long int location; // current location (position) in reference chromosome.
    char lastBase = 'X'; // the last base of reference line. this is for CG_only mode.
    UnsafeQueue<Position*> freePositionPool; // pool to store free position pointer for reference position.
    long long int refCoveredPosition; // this is the last position in reference chromosome we loaded in refPositions.
    ifstream refFile;
    ChromosomeFilePositions chromosomePos; // store the chromosome name and it's streamPos. To quickly find new chromosome in file.
    bool addedChrName = false;
    bool removedChrName = false;
		
		BufferedOutput out;

    Positions(string inputRefFileName, bool inputAddedChrName, bool inputRemovedChrName)
			: out(cout, 10000) {
        addedChrName = inputAddedChrName;
        removedChrName = inputRemovedChrName;
        refFile.open(inputRefFileName, ios_base::in);
        LoadChromosomeNamesPos();
    }

    ~Positions() {
        Position* pos;
        while(freePositionPool.popFront(pos)) {
            delete pos;
        }
    }

    /**
     * given the target Position output the corresponding position index in refPositions.
     */
    int getIndex(long long int &targetPos) {
        int firstPos = refPositions[0]->location;
        return targetPos - firstPos;
    }

    /**
     * given reference line (start with '>'), extract the chromosome information.
     * this is important when there is space in chromosome name. the SAM information only contain the first word.
     */
    string getChrName(string& inputLine) {
        string name;
        for (int i = 1; i < inputLine.size(); i++)
        {
            char c = inputLine[i];
            if (isspace(c)){
                break;
            }
            name += c;
        }

        if(removedChrName) {
            if(name.find("chr") == 0) {
                name = name.substr(3);
            }
        } else if(addedChrName) {
            if(name.find("chr") != 0) {
                name = string("chr") + name;
            }
        }
        return name;
    }


    /**
     * Scan the reference file. Record each chromosome and its position in file.
     */
    void LoadChromosomeNamesPos() {
        string line;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // this line is chromosome name
                chromosome = getChrName(line);
                streampos currentPos = refFile.tellg();
                chromosomePos.append(chromosome, currentPos);
            }
        }
        chromosomePos.sort();
        chromosome.clear();
    }

    /**
     * get a fasta line (not header), append the bases to positions.
     */
    void appendRefPosition(string& line) {
        
        // check the base one by one
        int len = line.size();

        refPositions.reserve(refPositions.size() + len);
        
        Position* newPos;

        #pragma unroll
        for (int i = 0; i < len; i++) {
            getFreePosition(newPos);
            newPos->set(curChromosomeId, location+i);
            char b = line[i];
            if (CG_only) {
                if (lastBase == 'C' && b == 'G') {
                    refPositions.back()->set('+');
                    newPos->set('-');
                }
            } else {
                if (b == convertFrom) {
                    newPos->set('+');
                } else if (b == convertFromComplement) {
                    newPos->set('-');
                }
            }
            refPositions.push_back(newPos);
            lastBase = b;
        }
        location += len;
    }

		void output_pos(const Position *pos)
		{
				if (!pos->isEmpty() && pos->strand != '?') {
						out.lock() << chromosomePos.getChromesomeString(pos->chromosomeId) << '\t'
											 << pos->location << '\t'
											 << pos->strand << '\t' 
											 << pos->convertedCount << '\t' 
											 << pos->unconvertedCount << '\n';
				}
		}

    /**
     * move the position which position smaller than refCoveredPosition - loadingBlockSize, output it.
     */
    void moveBlockToOutput() {
        if (refPositions.empty()) {
            return;
        }
        int index;
        int len = refPositions.size();
        for (index = 0; index < len; index++) {
            if (refPositions[index]->location < refCoveredPosition - loadingBlockSize) {
								output_pos(refPositions[index]);
								returnPosition(refPositions[index]);
            } else {
                break;
            }
        }
        if (index != 0) {
            refPositions.erase(refPositions.begin(), refPositions.begin()+index);
        }
    }

    /**
     * move all the refPosition into output pool.
     */
    void moveAllToOutput() {
        if (refPositions.empty()) {
            return;
        }
        for (int index = 0; index < refPositions.size(); index++) {
						output_pos(refPositions[index]);
						returnPosition(refPositions[index]);
        }
        refPositions.clear();
    }

    /**
     * initially load reference sequence for 2 million bp
     */
    void loadNewChromosome(string targetChromosome) {
        refFile.clear();
        // find the start position in file based on chromosome name.
        streampos startPos = chromosomePos.getChromosomePosInRefFile(targetChromosome);
        chromosome = targetChromosome;
        curChromosomeId = chromosomePos.findChromosome(targetChromosome, 0, chromosomePos.pos.size()-1);
        refFile.seekg(startPos, ios::beg);
        refCoveredPosition = 2 * loadingBlockSize;
        string line;
        lastBase = 'X';
        location = 0;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // this line is chromosome name
                return; // meet next chromosome, return it.
            } else {
                if (line.empty()) { continue; }
                // change all base to upper case
                for (int i = 0; i < line.size(); i++) {
                    line[i] = toupper(line[i]);
                }
                appendRefPosition(line);
                if (location >= refCoveredPosition) {
                    return;
                }
            }
        }
    }

    /**
     * load more Position (loadingBlockSize bp) to positions
     * if we meet next chromosome, return false. Else, return ture.
     */
    void loadMore() {
        refCoveredPosition += loadingBlockSize;
        string line;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // meet next chromosome, return.
                return ;
            } else {
                if (line.empty()) { continue; }

                // change all base to upper case
                for (int i = 0; i < line.size(); i++) {
                    line[i] = toupper(line[i]);
                }

                appendRefPosition(line);
                if (location >= refCoveredPosition) {
                    return ;
                }
            }
        }
    }


    /**
     * add position information from Alignment into ref position.
     */
    void appendPositions(Alignment& newAlignment) {
        if (!newAlignment.mapped || newAlignment.bases.empty()) {
            return;
        }
        long long int startPos = newAlignment.location; // 1-based position
        // find the first reference position in pool.
        int index = getIndex(newAlignment.location);

        for (int i = 0; i < newAlignment.sequence.size(); i++) {
            PosQuality* b = &newAlignment.bases[i];
            if (b->remove) {
                continue;
            }

            Position* pos = refPositions[index+b->refPos];
            assert (pos->location == startPos + b->refPos);

            if (pos->strand == '?') {
                // this is for CG-only mode. read has a 'C' or 'G' but not 'CG'.
                continue;
            }
            pos->appendBase(newAlignment.bases[i], newAlignment);
        }
    }

    /**
     * get a Position pointer from freePositionPool, if freePositionPool is empty, make a new Position pointer.
     */
    // 一次取一行
    void getFreePosition(Position*& newPosition) {
        if (freePositionPool.popFront(newPosition)) {
            return;
        } else {
            newPosition = new Position();
        }
    }

    /**
     * return the position to freePositionPool.
     */
    void returnPosition(Position* pos) {
        pos->initialize();
        freePositionPool.push(pos);
    }
};

#endif //POSITION_3N_TABLE_H
