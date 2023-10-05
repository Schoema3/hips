
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

class Integrator {
public:
    virtual void integrate(double dt) = 0;
    // Add other common integration methods here
    virtual ~Integrator() {} // Define a virtual destructor for proper cleanup
};

#endif

