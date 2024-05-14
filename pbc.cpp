#include <iostream>
#include <cmath>

#define L 5

int main()
{
    float x = -10.0;
    while (x < 10.0)
    {
        std::cout << x << " " << x - rint(x/L)*L << std::endl;
        x += 1.0;
    }
}