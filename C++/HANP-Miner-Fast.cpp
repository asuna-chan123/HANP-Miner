#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstring>
#include <algorithm>
#include <windows.h>

using namespace std;

#define K 6000     // Maximum number of sequences
#define M 1000     // Maximum pattern length
#define N 200000   // Maximum sequence length

struct sub_ptn_struct
{
    char start, end;
    int min, max;
};

struct seqdb
{
    int id;
    int len;
    char S[N];
} sDB[K];

int NumbS = 0;
int mingap, maxgap;
double minsup;

double unum[K];
int ww = 0;

sub_ptn_struct sub_ptn[M];
int ptn_len = 0;

vector<string> freArr[M];
vector<string> canArr[M];
vector<string> candidate;
int compnum = 0, frenum = 0;

double utility_val[256];

// 
int last_pos[M];

void init_utility()
{
    for (int i = 0; i < 256; i++) utility_val[i] = 0;
    utility_val['a'] = 3;
    utility_val['g'] = 6;
    utility_val['c'] = 6;
    utility_val['t'] = 4;
}

void deal_range(const string &pattern)
{
    ptn_len = 0;
    int len = (int)pattern.size();
    if (len == 1)
    {
        sub_ptn[0].start = pattern[0];
        sub_ptn[0].min = sub_ptn[0].max = 0;
        return;
    }
    for (int i = 0; i < len - 1; i++)
    {
        sub_ptn[ptn_len].start = pattern[i];
        sub_ptn[ptn_len].end = pattern[i + 1];
        sub_ptn[ptn_len].min = mingap;
        sub_ptn[ptn_len].max = maxgap;
        ptn_len++;
    }
}

int binary_search(int level, const string &cand, int low, int high)
{
    if (low > high) return -1;
    int prefix_len = level - 1;
    while (low <= high)
    {
        int mid = (low + high) / 2;
        int result = cand.compare(0, prefix_len, canArr[level - 1][mid], 0, prefix_len);
        if (result == 0)
        {
            int slow = low;
            int shigh = mid;
            if (cand.compare(0, prefix_len, canArr[level - 1][low], 0, prefix_len) == 0)
            {
                return low;
            }
            while (slow < shigh)
            {
                int smid = (slow + shigh) / 2;
                if (cand.compare(0, prefix_len, canArr[level - 1][smid], 0, prefix_len) == 0)
                {
                    shigh = smid;
                }
                else
                {
                    slow = smid + 1;
                }
            }
            return slow;
        }
        else if (result < 0)
        {
            high = mid - 1;
        }
        else
        {
            low = mid + 1;
        }
    }
    return -1;
}

void gen_candidate(int level)
{
    int size = (int)canArr[level - 1].size();
    int start = 0;
    int sub_len = level - 1;

    for (int i = 0; i < size; i++)
    {
        string R = canArr[level - 1][i].substr(1, sub_len);
        string Q = canArr[level - 1][start].substr(0, sub_len);
        if (Q != R)
        {
            start = binary_search(level, R, 0, size - 1);
        }
        if (start < 0 || start >= size)
        {
            start = 0;
        }
        else
        {
            while (start < size && canArr[level - 1][start].compare(0, sub_len, R, 0, sub_len) == 0)
            {
                string cand = canArr[level - 1][i];
                cand += canArr[level - 1][start][level - 1];
                candidate.push_back(cand);
                start++;
            }
            if (start >= size) start = 0;
        }
    }
}

void min_freItem()
{
    int counter[256] = {0};
    for (int t = 0; t < NumbS; t++)
    {
        const char *s = sDB[t].S;
        int len = sDB[t].len;
        for (int i = 0; i < len; i++)
        {
            counter[(unsigned char)s[i]]++;
        }
    }

    for (int c = 0; c < 256; c++)
    {
        if (counter[c] > 0)
        {
            string item(1, (char)c);
            double hupval = utility_val[c];
            if (counter[c] * hupval < minsup)
            {
                double uphupval = counter[c] * 6.0;
                if (uphupval >= minsup)
                {
                    canArr[0].push_back(item);
                }
            }
            else
            {
                freArr[0].push_back(item);
                canArr[0].push_back(item);
                unum[ww++] = hupval;
            }
        }
    }
}

int create_subnettree(const char *S, int len_S, int parent, int L)
{
    if (L > ptn_len + 1)
    {
        return 1;
    }
    int left = parent + sub_ptn[L - 2].min + 1;
    int right = min(parent + sub_ptn[L - 2].max + 1, len_S - 1);
    char target = sub_ptn[L - 2].end;

    for (int j = left; j <= right; j++)
    {
        if (S[j] == target)
        {
            if (j > last_pos[L - 1])
            {
                last_pos[L - 1] = j;
                if (create_subnettree(S, len_S, j, L + 1) == 1)
                {
                    return 1;
                }
            }
        }
    }
    return 0;
}

int create_nettree(const char *S, int len_S)
{
    int occurnum = 0;
    for (int i = 0; i <= ptn_len; i++)
    {
        last_pos[i] = -1;
    }

    int limit = len_S - ptn_len;
    char first_ch = sub_ptn[0].start;

    for (int i = 0; i < limit; i++)
    {
        if (S[i] != first_ch) continue;
        last_pos[0] = i;
        if (create_subnettree(S, len_S, i, 2) == 1)
        {
            occurnum++;
        }
    }
    return occurnum;
}

void read_file(const char *namefile)
{
    ifstream file(namefile);
    if (!file.is_open())
    {
        cerr << "Cannot open file: " << namefile << endl;
        exit(1);
    }
    string buff;
    int i = 0;
    while (getline(file, buff))
    {
        if (buff.empty()) continue;
        strcpy(sDB[i].S, buff.c_str());
        sDB[i].len = (int)buff.size();
        sDB[i].id = i + 1;
        i++;
    }
    NumbS = i;
}

int main(int argc, const char *argv[])
{
    if (argc < 5)
    {
        cout << "Usage: " << argv[0] << " <file> <mingap> <maxgap> <minsup>" << endl;
        return 1;
    }

    read_file(argv[1]);
    mingap = atoi(argv[2]);
    maxgap = atoi(argv[3]);
    minsup = atof(argv[4]);
    init_utility();

    int cannum = 0;
    DWORD begintime = GetTickCount();

    min_freItem();
    int f_level = 1;
    gen_candidate(f_level);

    while (!candidate.empty())
    {
        int cand_size = (int)candidate.size();
        for (int i = 0; i < cand_size; i++)
        {
            cannum++;
            compnum++;
            const string &p = candidate[i];

            double p_utility_sum = 0;
            for (char c : p)
            {
                p_utility_sum += utility_val[(unsigned char)c];
            }

            deal_range(p);

            int occnum = 0;
            for (int t = 0; t < NumbS; t++)
            {
                if (sDB[t].len >= ptn_len + 1)
                {
                    occnum += create_nettree(sDB[t].S, sDB[t].len);
                }
            }

            double hupval = (double)occnum * p_utility_sum / p.size();
            if (hupval >= minsup)
            {
                freArr[f_level].push_back(p);
                canArr[f_level].push_back(p);
                unum[ww++] = hupval;
            }
            else
            {
                double uphupval = (double)occnum * 6.0;
                if (uphupval >= minsup)
                {
                    canArr[f_level].push_back(p);
                }
            }
        }

        f_level++;
        candidate.clear();
        gen_candidate(f_level);
    }

    DWORD endtime = GetTickCount();

    for (int i = 0; i < f_level; i++)
    {
        for (size_t j = 0; j < freArr[i].size(); j++)
        {
            cout << freArr[i][j] << "   ";
            frenum++;
        }
        cout << endl;
    }

    cout << "The number of frequent patterns:" << frenum << endl;
    cout << "The time-consuming:" << endtime - begintime << "ms. \n";
    cout << "The number of calculation:" << compnum << " \n";

    time_t t = time(0);
    char tmp[32] = {0};
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S", localtime(&t));

    ofstream mcfile("result_file.txt", ios::app);
    mcfile << tmp << "\t" << argv[1] << "\t" << minsup << "\t\t[" << mingap << "," << maxgap << "]\t"
           << endtime - begintime << "\t\t\t" << frenum << "\t\t\t" << cannum << "\n";
    mcfile.close();

    ofstream maxfrefile("result_fature.txt", ios::app);
    maxfrefile << tmp << "\t" << argv[1] << "\t" << minsup << "\t\t[" << mingap << "," << maxgap << "]\t"
               << endtime - begintime << "\t\t\t" << frenum << "\n";

    for (int m = 0; m < f_level; m++)
    {
        for (size_t j = 0; j < freArr[m].size(); j++)
        {
            maxfrefile << freArr[m][j] << "\t";
        }
    }
    maxfrefile << "\n";
    for (int q = 0; q < ww; q++)
    {
        maxfrefile << unum[q] << "\t";
    }
    maxfrefile << "\n";
    maxfrefile.close();

    return 0;
}
