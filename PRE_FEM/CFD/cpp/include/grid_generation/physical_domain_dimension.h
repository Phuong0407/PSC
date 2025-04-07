#ifndef PHYSICAL_DOMAIN_DIMENSION_H
#define PHYSICAL_DOMAIN_DIMENSION_H

class PhysicalDomainDimension{
private:
    void initVerticalDimension(double L1, double L2);
    void initHorizontalDimension(double H1, double H2, double H3);
    
public:
    double L1, L2;
    double H1, H2, H3;
    PhysicalDomainDimension() = default;
    void initPhysicalDomainDimension(double L1, double L2, double H1, double H2, double H3);
};

#endif