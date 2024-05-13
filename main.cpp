#include <iostream>
#include <random>
#include <cmath>
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

#define T 300.0
#define KB 8.31e-3

struct float3
{
    float x, y, z;
};

struct Atom
{
    int index;
    std::string name;

    float q;
    float m;
    float sigma;
    float eps;
};

void saveCoordinates(const std::string filename,
                   const std::vector<Atom>& atoms,
                   const std::vector<float3>& r,
                   const std::vector<float3>& v)
{
    FILE* fout = fopen(filename.c_str(), "w");
    fprintf(fout, "NaCl\n%d\n", N);
    for (int i = 0; i < N; i++)
    {
        fprintf(fout, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
            i, "NaCl",
            atoms[i].name.c_str(), atoms[i].index,
            r[i].x, r[i].y, r[i].z,
            v[i].x, v[i].y, v[i].z
        );
    }
    fprintf(fout, "%8.3f %8.3f %8.3f\n\n", L, L, L);
    fclose(fout);
}

int main()
{
    std::vector<Atom> atoms(N);

    std::vector<float3> r(N);
    std::vector<float3> v(N);
    std::vector<float3> f(N);
    
    for (int i = 0; i < N; i++)
    {
        if (i < N/2)
        {
            atoms[i].index = i;
            atoms[i].name = "Na";
            atoms[i].q = Q_NA;
            atoms[i].m = M_NA;
            atoms[i].sigma = SIGMA_NA;
            atoms[i].eps = EPSILON_NA;
        }
        else
        {
            atoms[i].index = i;
            atoms[i].name = "Cl";
            atoms[i].q = Q_CL;
            atoms[i].m = M_CL;
            atoms[i].sigma = SIGMA_CL;
            atoms[i].eps = EPSILON_CL;
        }
    }

    std::mt19937 randomGenerator(123456);
    std::uniform_real_distribution<> distributionX(0.0, L);
    std::normal_distribution<> distributionV(0.0, sqrtf(KB*T));

    for (int i = 0; i < N; i++)
    {
        r[i].x = distributionX(randomGenerator);
        r[i].y = distributionX(randomGenerator);
        r[i].z = distributionX(randomGenerator);

        v[i].x = distributionV(randomGenerator)/sqrtf(atoms[i].m);
        v[i].y = distributionV(randomGenerator)/sqrtf(atoms[i].m);
        v[i].z = distributionV(randomGenerator)/sqrtf(atoms[i].m);
    }

    saveCoordinates("coord.gro", atoms, r, v);
}

