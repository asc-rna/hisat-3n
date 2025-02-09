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
#include <shared_mutex>

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
    // mutex mutex_;
public:
    // string chromosome; // reference chromosome name
    long long int location; // 1-based position
    short chromosomeId;
    unsigned short convertedCount = 0;
    unsigned short unconvertedCount = 0;
    
    char strand; // +(REF) or -(REF-RC)
    // string convertedQualities; // each char is a mapping quality on this position for converted base.
    // string unconvertedQualities; // each char is a mapping quality on this position for unconverted base.

    bool empty = true;
    // vector<uniqueID> uniqueIDs; // each value represent a readName which contributed the base information.
    //                           // readNameIDs is to make sure no read contribute 2 times in same position.

    void initialize() {
        // chromosome.clear();
        chromosomeId = -1;
        location = -1;
        strand = '?';
        // convertedQualities.clear();
        // unconvertedQualities.clear();
        convertedCount = 0;
        unconvertedCount = 0;
        empty = true;
        // vector<uniqueID>().swap(uniqueIDs);
    }

    Position(){
        initialize();
    };

    /**
     * return true if there is mapping information in this reference position.
     */
    inline bool isEmpty() {
        // return convertedQualities.empty() && unconvertedQualities.empty();
        return empty;
    }

    /**
     * set the chromosome, location (position), and strand information.
     */

    inline void set (int InChromosomeId, long long int inputLoc) {
        chromosomeId = InChromosomeId;
        location = inputLoc + 1;
    }

    inline void set(char inputStrand) {
        strand = inputStrand;
    }

    /**
     * binary search of readNameID in readNameIDs.
     * always return a index.
     * if cannot find, return the index which has bigger value than input readNameID.
     */
    // int searchReadNameID (unsigned long long&readNameID, int start, int end) {
    //     if (uniqueIDs.empty()) {
    //         return 0;
    //     }
    //     if (start <= end) {
    //         int middle = (start + end) / 2;
    //         if (uniqueIDs[middle].readNameID == readNameID) {
    //             return middle;
    //         }
    //         if (uniqueIDs[middle].readNameID > readNameID) {
    //             return searchReadNameID(readNameID, start, middle-1);
    //         }
    //         return searchReadNameID(readNameID, middle+1, end);
    //     }
    //     return start; // return the bigger one
    // }


    /**
     * with a input readNameID, add it into readNameIDs.
     * if the input readNameID already exist in readNameIDs, return false.
     */
    // bool appendReadNameID(PosQuality& InBase, Alignment& InAlignment) {
    //     int idCount = uniqueIDs.size();
    //     if (idCount == 0 || InAlignment.readNameID > uniqueIDs.back().readNameID) {
    //         uniqueIDs.emplace_back(InAlignment.readNameID, InBase.converted, InBase.qual);
    //         return true;
    //     }
    //     int index = searchReadNameID(InAlignment.readNameID, 0, idCount);
    //     if (uniqueIDs[index].readNameID == InAlignment.readNameID) {
    //         // if the new base is consistent with exist base's conversion status, ignore
    //         // otherwise, delete the exist conversion status
    //         if (uniqueIDs[index].removed) {
    //             return false;
    //         }
    //         if (uniqueIDs[index].isConverted != InBase.converted) {
    //             uniqueIDs[index].removed = true;
    //             if (uniqueIDs[index].isConverted) {
    //                 for (int i = 0; i < convertedQualities.size(); i++) {
    //                     if (convertedQualities[i] == InBase.qual) {
    //                         convertedQualities.erase(convertedQualities.begin()+i);
    //                         return false;
    //                     }
    //                 }
    //             } else {
    //                 for (int i = 0; i < unconvertedQualities.size(); i++) {
    //                     if (unconvertedQualities[i] == InBase.qual) {
    //                         unconvertedQualities.erase(unconvertedQualities.begin()+i);
    //                         return false;
    //                     }
    //                 }
    //             }
    //         }
    //         return false;
    //     } else {
    //         uniqueIDs.emplace(uniqueIDs.begin()+index, InAlignment.readNameID, InBase.converted, InBase.qual);
    //         return true;
    //     }
    // }

    /**
     * append the SAM information into this position.
     */
    void appendBase (PosQuality& input, Alignment& a) {
        // mutex_.lock();
            if (empty) empty = false;
            if (input.converted) {
                // convertedQualities += input.qual;
                convertedCount++;
            } else {
                // unconvertedQualities += input.qual;
                unconvertedCount++;
            }
        
        // mutex_.unlock();
    }
};

/**
 * store all reference position in this class.
 */
class Positions{
public:
    // vector<Position*> refPositions; // the pool of all current reference position.
    vector<Position> ToOutRefPositions; // the pool of all current reference position.
    vector<Position> refPositions;

    string chromosome; // current reference chromosome name.'
    int curChromosomeId;
    long long int location; // current location (position) in reference chromosome.
    char lastBase = 'X'; // the last base of reference line. this is for CG_only mode.
    SafeQueue<string*> linePool; // pool to store unprocessed SAM line.
    SafeQueue<string*> freeLinePool; // pool to store free string pointer for SAM line.
    UnsafeQueue<Position*> freePositionPool; // pool to store free position pointer for reference position.
    SafeQueue<Position*> outputPositionPool; // pool to store the reference position which is loaded and ready to output.
    bool working;
    mutex mutex_;
    mutable shared_mutex s_mutex;
    long long int refCoveredPosition; // this is the last position in reference chromosome we loaded in refPositions.
    ifstream refFile;
    vector<mutex*> workerLock; // one lock for one worker thread.
    int nThreads = 1;
    ChromosomeFilePositions chromosomePos; // store the chromosome name and it's streamPos. To quickly find new chromosome in file.
    bool addedChrName = false;
    bool removedChrName = false;


    bool output = false;
    bool finalOut = false;

    thread* outputThread;
    ofstream tableFile;

    Alignment tmpAlignment;

    
    
    

    Positions(string inputRefFileName, int inputNThreads, bool inputAddedChrName, bool inputRemovedChrName,string outputFileName) {
        working = true;
        nThreads = inputNThreads;
        addedChrName = inputAddedChrName;
        removedChrName = inputRemovedChrName;
        for (int i = 0; i < nThreads; i++) {
            workerLock.push_back(new mutex);
        }
        refFile.open(inputRefFileName, ios_base::in);
        LoadChromosomeNamesPos();

        if (!outputFileName.empty()) {
            tableFile.open(outputFileName, ios_base::out );
            tableFile << "ref\tpos\tstrand\tconvertedBaseCount\tunconvertedBaseCount\n";
        }
    }

    ~Positions() {
        for (int i = 0; i < workerLock.size(); i++) {
            delete workerLock[i];
        }

        tableFile.close();

        Position* pos;
        while(freePositionPool.popFront(pos)) {
            delete pos;
        }
    }

    void startOutput(bool final_ = false){



        if (refPositions.empty()) {
            return;
        }


        // 写锁保护
        {

        std::lock_guard<std::shared_mutex> lock(s_mutex);
        if (!ToOutRefPositions.empty()) {

        for (Position& pos:ToOutRefPositions){
            pos.initialize();
        }
        auto start = ToOutRefPositions.begin() ;
        auto end =  ToOutRefPositions.end();
        refPositions.insert(refPositions.end(), std::make_move_iterator(start), std::make_move_iterator(end));    
        ToOutRefPositions.erase(start, end);
        }
        

        if (!final_){

        auto start = refPositions.begin() ;
        auto end =  refPositions.begin() + loadingBlockSize;
        // 使用 std::move 移动元素到 destination
        ToOutRefPositions.insert(ToOutRefPositions.end(), std::make_move_iterator(start), std::make_move_iterator(end));
        // 清除原 vector 中的这部分元素
        refPositions.erase(start, end);

        // printf("ToOutRefPositions size: %d\n", ToOutRefPositions.size());
        // printf("refPositions size: %d\n", refPositions.size());

        // for (int i=0;i<ToOutRefPositions.size();i++){
        //     printf("ToOutRefPositions[%d]: %d\n", i, ToOutRefPositions[i]->location);
        //     delete ToOutRefPositions[i];
        // }
        // ToOutRefPositions.clear();
        } else {
            std::swap(refPositions, ToOutRefPositions);
            refPositions.clear();
        }

        }

        long outLimitPos = refCoveredPosition - loadingBlockSize;


        outputOnce(outLimitPos, final_);

        if (final_) ToOutRefPositions.clear();


    }


    /**
     * given the target Position output the corresponding position index in refPositions.
     */
    int getIndex(long long int &targetPos) {
        int firstPos = refPositions[0].location;
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
        // cout << "chromosome name: " << name << endl;
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
    inline void appendRefPosition(string& line,int& cur) {
        
        // check the base one by one
        int len = line.size();


        Position* newPos;

        #pragma unroll(60)
        for (int i = 0; i < len; i++) {
            
            refPositions[cur + i].set(curChromosomeId, location+i);
            char b = line[i];
            // if (CG_only) {
            //     if (lastBase == 'C' && b == 'G') {
            //         refPositions[.back()]->set('+');
            //         newPos->set('-');
            //     }
            // } else {
                if (b == convertFrom) {
                    refPositions[cur+i].set('+');
                } else if (b == convertFromComplement) {
                    refPositions[cur+i].set('-');
                }
            // }
            lastBase = b;
        }
        location += len;
        cur += len;
    }

    /**
     * if we can go through all the workerLock, that means no worker is appending new position.
     */
    void appendingFinished() {
        for (int i = 0; i < nThreads; i++) {
            workerLock[i]->lock();
            workerLock[i]->unlock();
        }
    }

    /**
     * the output function for output thread.
     */

    inline void outputPosition(Position* pos,ostream* out_) {

    }

    int outputOnce(long outLimitPos, bool final_=false) {
        if (ToOutRefPositions.empty()) {
            return 0;
        }
        int i;
        for (i = 0; i < ToOutRefPositions.size(); i++){
            Position& pos = ToOutRefPositions[i];

            if (final_ || pos.location <= outLimitPos){

        
            if (pos.isEmpty() || pos.strand == '?') {

            }
            
            else {
            const string& chr = chromosomePos.getChromesomeString(pos.chromosomeId);
            tableFile << chr << '\t'
                << to_string(pos.location) << '\t'
                << pos.strand << '\t'
                << pos.convertedCount << '\t'
                << pos.unconvertedCount << '\n';


            }

            } else {
                break;
            }
        }

            return i;
    }


    // void outputFunction(string outputFileName) {
    //     ostream* out_ = &cout;
    //     out_ = &cout;
    //     ofstream tableFile;
    //     if (!outputFileName.empty()) {
    //         tableFile.open(outputFileName, ios_base::out );
    //         out_ = &tableFile;
    //     }

    //     // *out_ << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunconvertedBaseQualities\tunconvertedBaseCount\n";

    //     *out_ << "ref\tpos\tstrand\tconvertedBaseCount\tunconvertedBaseCount\n";
    //     while (working) {


    //         // if (outputPositionPool.popFront(pos)) {
    //         //     const string& chr = chromosomePos.getChromesomeString(pos->chromosomeId);
    //         //     *out_ << chr << '\t'
    //         //               << to_string(pos->location) << '\t'
    //         //               << pos->strand << '\t'
    //         //               << pos->convertedCount << '\t'
    //         //               << pos->unconvertedCount << '\n';
    //         //     delete pos;
    //         // } else {
    //         //     this_thread::sleep_for (std::chrono::microseconds(1));
    //         // }

    //         if (output){
    //             int i;
    //             for (i = 0; i < ToOutRefPositions.size(); i++){
    //                 Position* pos = ToOutRefPositions[i];

    //                 if (finalOut || pos->location < outLimitPos){

                
    //                 if (pos->isEmpty() || pos->strand == '?') {
    //                     delete pos;
    //                 }
                    
    //                 else {
    //                 const string& chr = chromosomePos.getChromesomeString(pos->chromosomeId);
    //                 *out_ << chr << '\t'
    //                     << to_string(pos->location) << '\t'
    //                     << pos->strand << '\t'
    //                     << pos->convertedCount << '\t'
    //                     << pos->unconvertedCount << '\n';
    //                 delete pos;
    //                 }

    //                 } else {
    //                     break;
    //                 }
    //             }
    //                     if (i != 0) {
    //                     // printf("i: %d\n", i);
    //                     ToOutRefPositions.erase(ToOutRefPositions.begin(), ToOutRefPositions.begin()+i);
    //     }

                
    //             finishOutput();
    //         }
    //         else {
    //             this_thread::sleep_for (std::chrono::microseconds(1));
    //         }
    //     }
    //     tableFile.close();
    // }

    /**
     * move the position which position smaller than refCoveredPosition - loadingBlockSize, output it.
     */
    // void moveBlockToOutput() {
    //     if (refPositions.empty()) {
    //         return;
    //     }
    //     int index;
    //     int len = refPositions.size();
    //     for (index = 0; index < len; index++) {
    //         if (refPositions[index].location < refCoveredPosition - loadingBlockSize) {
    //             if (refPositions[index].isEmpty() || refPositions[index].strand == '?') {
    //                 returnPosition(refPositions[index]);
    //             } else {
    //                 outputPositionPool.push(refPositions[index]);
    //             }
    //         } else {
    //             break;
    //         }
    //     }
    //     if (index != 0) {
    //         // printf("index: %d\n", index);
    //         refPositions.erase(refPositions.begin(), refPositions.begin()+index);

    //     }
    // }

    // /**
    //  * move all the refPosition into output pool.
    //  */
    // void moveAllToOutput() {
    //     if (refPositions.empty()) {
    //         return;
    //     }
    //     for (int index = 0; index < refPositions.size(); index++) {
    //         if (refPositions[index]->isEmpty() || refPositions[index]->strand == '?') {
    //             returnPosition(refPositions[index]);
    //         } else {
    //             outputPositionPool.push(refPositions[index]);
    //         }
    //     }
    //     refPositions.clear();
    // }

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

        refPositions.resize(3*loadingBlockSize);
        string line;
        lastBase = 'X';
        location = 0;
        int cur = 0;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // this line is chromosome name
                return; // meet next chromosome, return it.
            } else {
                if (line.empty()) { continue; }
                // change all base to upper case

                // for (int i = 0; i < line.size(); i++) {
                //     line[i] = toupper(line[i]);
                // }
                appendRefPosition(line,cur);
                if (location >= refCoveredPosition) {
                    // printf("cur %d\n", cur);
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
        int cur = loadingBlockSize;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // meet next chromosome, return.
                // printf("meet next chromosome\n");
                return ;
            } else {
                if (line.empty()) { continue; }

                // // change all base to upper case
                // for (int i = 0; i < line.size(); i++) {
                //     line[i] = toupper(line[i]);
                // }

                appendRefPosition(line,cur);
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

        int len = refPositions.size();

        for (int i = 0; i < newAlignment.sequence.size(); i++) {
            PosQuality* b = &newAlignment.bases[i];
            if (b->remove) {
                continue;
            }

            Position& pos = refPositions[index+b->refPos];
            assert (pos.location == startPos + b->refPos);

            if (pos.strand == '?') {
                // this is for CG-only mode. read has a 'C' or 'G' but not 'CG'.
                continue;
            }
            pos.appendBase(newAlignment.bases[i], newAlignment);
        }
    }

    /**
     * get a string pointer from freeLinePool, if freeLinePool is empty, make a new string pointer.
     */
    void getFreeStringPointer(string*& newLine) {
        if (freeLinePool.popFront(newLine)) {
            return;
        } else {
            newLine = new string();
        }
    }

    /**
     * get a Position pointer from freePositionPool, if freePositionPool is empty, make a new Position pointer.
     */
    // 一次取一行
    // void getFreePosition(Position*& newPosition) {
    //     // while (outputPositionPool.size() >= 10000) {
    //     //     this_thread::sleep_for (std::chrono::microseconds(1));
    //     // }
    //     if (freePositionPool.popFront(newPosition)) {
    //         return;
    //     } else {
    //         newPosition = new Position();
    //     }
    // }

    /**
     * return the line to freeLinePool
     */
    void returnLine(string* line) {
        line->clear();
        freeLinePool.push(line);
    }

    /**
     * return the position to freePositionPool.
     */
    void returnPosition(Position* pos) {
        pos->initialize();
        freePositionPool.push(pos);
    }

    /**
     * this is the working function.
     * it take the SAM line from linePool, parse it.
     */
    void append(int threadID) {
        string* line;
        Alignment newAlignment;
       

        while (working) {
             // 读锁保护
            std::shared_lock<std::shared_mutex> lock(s_mutex);
            workerLock[threadID]->lock();
            if(!linePool.popFront(line)) {
                workerLock[threadID]->unlock();
                this_thread::sleep_for (std::chrono::nanoseconds(1));
                continue;
            }
            while (refPositions.empty()) {
                this_thread::sleep_for (std::chrono::microseconds(1));
            }
            newAlignment.parse(line);
            returnLine(line);
            appendPositions(newAlignment);
            workerLock[threadID]->unlock();
        }
    }

    void appendSync(string* line) {
        tmpAlignment.parse(line);
        returnLine(line);
        appendPositions(tmpAlignment);

    }
};

#endif //POSITION_3N_TABLE_H
