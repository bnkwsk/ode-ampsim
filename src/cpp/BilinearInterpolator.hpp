#include <cmath>

// z = f(x, y); x - linear, y - linear
class BilinearInterpolator
{
    typedef double ScalarType;

    ScalarType (*function)(ScalarType, ScalarType);
    ScalarType **data;
    ScalarType minX, maxX, minY, maxY, stepX, stepY;
    int nX, nY;

    public:
        BilinearInterpolator(ScalarType (*_function)(ScalarType, ScalarType), ScalarType _minX, ScalarType _maxX, ScalarType _minY, ScalarType _maxY, ScalarType _stepX, ScalarType _stepY) : function(_function), minX(_minX), maxX(_maxX), minY(_minY), maxY(_maxY), stepX(_stepX), stepY(_stepY)
        {
            nX = ceil((maxX - minX) / stepX);
            nY = ceil((maxY - minY) / stepY);
            data = new ScalarType*[nX];
            int iX = 0;
            int iY = 0;
            for(ScalarType x = minX; x <= maxX; x += stepX)
            {
                data[iX] = new ScalarType[nY];
                iY = 0;
                for(ScalarType y = minY; y <= maxY; y += stepY)
                {
                    data[iX][iY] = function(x, y);
                    ++iY;
                }
                ++iX;
            }
        }

        ~BilinearInterpolator()
        {
            for(int i = 0; i < nX; ++i)
                delete[] data[i];
            delete[] data;
        }

        ScalarType operator()(ScalarType x, ScalarType y)
        {
            if(x < minX || x > maxX || y < minY || y > maxY || x != x || y != y)
                return function(x, y);
            int xLo = floor((x - minX) / stepX);
            int xHi = ceil((x - minX) / stepX);
            int yLo = floor((y - minY) / stepY);
            int yHi = ceil((y - minY) / stepY);
            ScalarType tX = (x - minX) / stepX - xLo;
            ScalarType tY = (y - minY) / stepY -  yLo;

            return data[xLo][yLo] * (1.0 - tX) * (1.0 - tY) +
                   data[xHi][yLo] * tX * (1.0 - tY) +
                   data[xLo][yHi] * (1.0 - tX) * tY +
                   data[xHi][yHi] * tX * tY;
        }
};