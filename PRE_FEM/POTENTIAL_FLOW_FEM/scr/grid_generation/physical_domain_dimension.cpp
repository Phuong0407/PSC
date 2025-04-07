#include "../include/grid_generation/physical_domain_dimension.h"

void PhysicalDomainDimension::initVerticalDimension(double L1, double L2) {
    this->L1 = L1;
    this->L2 = L2;
}

void PhysicalDomainDimension::initHorizontalDimension(double H1, double H2, double H3) {
    this->H1 = H1;
    this->H2 = H2;
    this->H3 = H3;
}

void PhysicalDomainDimension::initPhysicalDomainDimension(double L1, double L2, double H1, double H2, double H3) {
    initVerticalDimension(L1, L2);
    initHorizontalDimension(H1, H2, H3);
}