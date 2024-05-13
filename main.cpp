#include <iostream>
#include <random>
#include <vector>

#define L 5.0
#define N 100

#define Q_NA 1.0
#define Q_CL -1.0
#define M_NA 22.98
#define M_CL 35.45

//   SOD    11    22.989770    0.000     A  0.251367073323  0.1962296  ; Sodium Ion
//   CLA    17    35.450000    0.000     A  0.404468018036  0.6276000  ; Chloride Ion

#define SIGMA_NA 0.2514
#define SIGMA_CL 0.4045
#define EPSILON_NA 0.1962
#define EPSILON_CL 0.6276

#define COULOMB 138.9118

struct float3
{
    float x, y, z;
};

void savePositions(const std::string filename, const std::vector<float3>& r)
{
    FILE* fout = fopen(filename.c_str(), "w");
    fprintf(fout, "NaCl\n%d\n", N);
    for (int i = 0; i < N; i++)
    {
        fprintf(fout, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
            i, "NaCl",
            "At", i,
            r[i].x, r[i].y, r[i].z,
            0.0, 0.0, 0.0
        );
    }
    fprintf(fout, "%8.3f %8.3f %8.3f\n\n", L, L, L);
    fclose(fout);
}

int main()
{
    std::vector<float3> r(N);
    std::mt19937 randomGenerator(123456);
    std::uniform_real_distribution<> distributionX(0, L);

    for (int i = 0; i < N; i++)
    {
        r[i].x = distributionX(randomGenerator);
        r[i].y = distributionX(randomGenerator);
        r[i].z = distributionX(randomGenerator);
    }

    savePositions("pos.gro", r);
}

