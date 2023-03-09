#ifndef _LACT_RUNPARA_H
#define _LACT_RUNPARA_H

#include <string>
#include <vector>
class LACT_RUNPARA
{
    public:
        std::string out_file;
        std::string lookup_file;
        bool havelookup;

        bool writelookup;

        bool DrawMode;

        int special_entries;
        bool SelectTel;
        bool ResetWeight;
        double num_weight[4]{1};
        double max_dist;
        double min_tel;
        std::vector<int> Only_Telescope;
        std::vector<int> Draw_Events;
        std::vector<std::string> input_file;

        LACT_RUNPARA();
        ~LACT_RUNPARA();

        void ProcessCommandLine(int argc, char** argv);
        bool ProcessCommandLine(int argc, char** argv, bool lookup_flag);
        std::string  GetLookupName()
        {
            return lookup_file;
        }
        int GetOnlyNum()
        {
            return Only_Telescope.size();
        }
        std::vector<int> GetOnlyTels()
        {
            return Only_Telescope;
        }
};


















#endif