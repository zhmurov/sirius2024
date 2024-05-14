#include <iostream>
#include <random>
#include <cmath>
#include <vector>

#define L 5.0
#define N 100

#define tau 0.001
#define relax 0.1
#define NSTEPS 10000
#define STRIDE 100

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
                     const std::string modifier,
                     const std::vector<Atom>& atoms,
                     const std::vector<float3>& r,
                     const std::vector<float3>& v)
{
    FILE* fout = fopen(filename.c_str(), modifier.c_str());
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
    fprintf(fout, "%8.3f %8.3f %8.3f\n", L, L, L);
    fclose(fout);
}

float transferPBC(float x)
{
    while (x < 0)
    {
        x += L;
    }
    while (x > L)
    {
        x -= L;
    }
    return x;
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

    std::mt19937 randomGenerator(34443);
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

        f[i].x = 0.0;
        f[i].y = 0.0;
        f[i].z = 0.0;
    }

    saveCoordinates("coord.gro", "w", atoms, r, v);

    double temperature = 0;
    float scaleV = 1.0;
    //Integrate equations of motion
    for(int step = 0; step <= NSTEPS; step++)
    {
        // Compute force
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i; j++)
            {
                float dx = r[i].x - r[j].x;
                float dy = r[i].y - r[j].y;
                float dz = r[i].z - r[j].z;

                dx -= rint(dx/L)*L;
                dy -= rint(dy/L)*L;
                dz -= rint(dz/L)*L;
                
                // Coulomb
                float dr2 = dx*dx + dy*dy + dz*dz;
                float dr = sqrtf(dr2);
                float dfc = COULOMB*atoms[i].q*atoms[j].q/(dr2*dr);

                // VdW
                float sigma = 0.5*(atoms[i].sigma + atoms[j].sigma);
                float eps = -12.0*sqrtf(atoms[i].eps*atoms[j].eps);
                float sor2 = (sigma*sigma)/dr2;
                float sor6 = sor2*sor2*sor2;
                float dflj = eps*sor6*(1.0 - sor6)/dr2;

                float df = dfc + dflj;

                f[i].x += df*dx;
                f[i].y += df*dy;
                f[i].z += df*dz;

                f[j].x += -df*dx;
                f[j].y += -df*dy;
                f[j].z += -df*dz;
            }
        }

        for (int i = 0; i < N; i++)
        {
            // m*dv/dt = f:
            v[i].x = v[i].x + tau*f[i].x/atoms[i].m;
            v[i].y = v[i].y + tau*f[i].y/atoms[i].m;
            v[i].z = v[i].z + tau*f[i].z/atoms[i].m;

            f[i].x = 0.0;
            f[i].y = 0.0;
            f[i].z = 0.0;

            // dx/dt = v:
            r[i].x = r[i].x + tau*v[i].x;
            r[i].y = r[i].y + tau*v[i].y;
            r[i].z = r[i].z + tau*v[i].z;

            r[i].x = transferPBC(r[i].x);
            r[i].y = transferPBC(r[i].y);
            r[i].z = transferPBC(r[i].z);

            temperature += atoms[i].m*
            (v[i].x*v[i].x + v[i].y*v[i].y + v[i].z*v[i].z);

            
            //Energy minimization
            if (step < 1000)
            {
                v[i].x = 0.0;
                v[i].y = 0.0;
                v[i].z = 0.0;
            } else
            if (step >= 1000 && step < 2000)
            {
                v[i].x = distributionV(randomGenerator)/sqrtf(atoms[i].m);
                v[i].y = distributionV(randomGenerator)/sqrtf(atoms[i].m);
                v[i].z = distributionV(randomGenerator)/sqrtf(atoms[i].m);
            } else {
                v[i].x *= scaleV;
                v[i].y *= scaleV;
                v[i].z *= scaleV;
            }
            
        }

        if (step % STRIDE == 0)
        {
            temperature /= (N*STRIDE*KB*3);
            std::cout << temperature << std::endl;
            scaleV = sqrtf(1.0 - ((temperature - T)/T)*(tau/relax));
            temperature = 0.0;
            saveCoordinates("coord.gro", "a", atoms, r, v);
        }
    }

    
}

