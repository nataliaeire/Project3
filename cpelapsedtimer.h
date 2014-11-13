#ifndef CPELAPSEDTIMER_H
#define CPELAPSEDTIMER_H
#include <time.h>

class CPTimingObject {
private:
    double m_timeElapsed;
    clock_t m_startedAt;
public:
    CPTimingObject() : m_timeElapsed(0), m_startedAt(0) { }

    void start() {
         m_startedAt = clock();
    }

    double stop() {
        double t = double(clock() - m_startedAt)/CLOCKS_PER_SEC;
        m_timeElapsed += t;
        return t;
    }

    double elapsedTime() { return m_timeElapsed; }
};

class CPElapsedTimer
{
public:
    CPElapsedTimer();

    static CPElapsedTimer& getInstance()
    {
        static CPElapsedTimer instance; // Guaranteed to be destroyed.
                                 // Instantiated on first use.
        return instance;
    }

    clock_t        m_startedAt;
    CPTimingObject m_calculateForces;
    CPTimingObject m_sortBodies;
    CPTimingObject m_gatherEnergies;

    static CPTimingObject &calculateForces() { return CPElapsedTimer::getInstance().m_calculateForces; }
    static CPTimingObject &sortBodies() { return CPElapsedTimer::getInstance().m_sortBodies; }
    static CPTimingObject &gatherEnergies() { return CPElapsedTimer::getInstance().m_gatherEnergies; }

    static double totalTime() { return double(clock() - CPElapsedTimer::getInstance().m_startedAt)/ CLOCKS_PER_SEC; }
};

#endif // CPELAPSEDTIMER_H
