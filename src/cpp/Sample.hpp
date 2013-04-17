#ifndef SAMPLE_HPP
#define SAMPLE_HPP

class Sample
{
    double value;
    bool last;

    public:
        Sample(double v=0.0, bool l=false) : value(v), last(l)
        {
            //
        }

        double getValue()
        {
            return value;
        }

        bool isLast()
        {
            return last;
        }
};

#endif