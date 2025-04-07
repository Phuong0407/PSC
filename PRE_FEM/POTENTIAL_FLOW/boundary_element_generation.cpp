#include "boundary_element_generation.h"
#include <algorithm>

void BoundaryElementGeneration::initVerticalPhysicalDimension(double L1, double L2) {
    this->L1 = L1;
    this->L2 = L2;
}
void BoundaryElementGeneration::initHorizontalPhysicalDimension(double H1, double H2, double H3) {
    this->H1 = H1;
    this->H2 = H2;
    this->H3 = H3;
}

void BoundaryElementGeneration::initPhysicalDomainDimension(double L1, double L2, double H1, double H2, double H3) {
    initVerticalPhysicalDimension(L1, L2);
    initHorizontalPhysicalDimension(H1, H2, H3);
}

void BoundaryElementGeneration::initBoundaryElementDimension(Size NumElemFirstBound, Size NumElemSecondBound, Size NumElemThirdBound, Size NumElemFourthBound, Size NumElemFifthBound, Size NumElemSixthBound) {
    this->NumElemFirstBound = NumElemFirstBound;
    this->NumElemSecondBound = NumElemSecondBound;
    this->NumElemThirdBound = NumElemThirdBound;
    this->NumElemFourthBound = NumElemFourthBound;
    this->NumElemFifthBound = NumElemFifthBound;
    this->NumElemSixthBound = NumElemSixthBound;
    this->NumElemTot = this->NumElemFirstBound + this->NumElemSecondBound + this->NumElemThirdBound + this->NumElemFourthBound + this->NumElemFifthBound + this->NumElemSixthBound;
    this->NumNodeTot = this->NumElemTot;
}

void BoundaryElementGeneration::buildElemConn() {
    ElemConn.clear();
    ElemConn.reserve(NumNodeTot);
    for (std::size_t i = 0; i < NumNodeTot; ++i) {
        ElemConn.emplace_back(i, (i + 1) % NumNodeTot);
    }
}

void BoundaryElementGeneration::buildNodeCoord() {
    double dz = H2 / NumElemFirstBound;
    double a = (L1 - L2)/H2, b = - a * H1;
    for (std::size_t i = 0; i < NumElemFirstBound; ++i) {
        double z = H1 + dz * i;
        double r = a * z + b;
        NodeCoord.emplace_back(z, r);
    }
    dz = H3 / NumElemSecondBound;
    for (std::size_t i = 0; i < NumElemSecondBound; ++i) {
        double z = H1 + H2 + dz * i;
        NodeCoord.emplace_back(z, L1 - L2);
    }
    double dr = L2 / NumElemThirdBound;
    for (std::size_t i = 0; i < NumElemThirdBound; ++i) {
        double r = L1 - L2 + dr * i;
        NodeCoord.emplace_back(H1 + H2 + H3, r);
    }
    dz = (H1 + H2 + H3) / NumElemFourthBound;
    for (std::size_t i = 0; i < NumElemFourthBound; ++i) {
        double z = H1 + H2 + H3 - dz * i;
        NodeCoord.emplace_back(z, L1);
    }
    dr = L1 / NumElemFifthBound;
    for (std::size_t i = 0; i < NumElemFifthBound; ++i) {
        double r = L1 - dr * i;
        NodeCoord.emplace_back(0, r);
    }
    dz = H1 / NumElemSixthBound;
    for(std::size_t i = 0; i < NumElemSixthBound; ++i) {
        double z = dz * i;
        NodeCoord.emplace_back(z, 0);
    }
}

void BoundaryElementGeneration::buildOutwardNormVec() {
    OutwardNormVec.clear();
    OutwardNormVec.reserve(NumElemTot);

    for (std::size_t i = 0; i < NumElemTot; ++i) {
        double z1 = NodeCoord[i].z, r1 = NodeCoord[i].r;
        double z2 = NodeCoord[(i + 1) % NumNodeTot].z, r2 = NodeCoord[(i + 1) % NumNodeTot].r;
        // double z2 = NodeCoord[i + 1].z, r2 = NodeCoord[i + 1].r;

        double dz = z2 - z1, dr = r2 - r1;
        double normZ = -dr, normR = dz;

        double length = std::sqrt(normZ * normZ + normR * normR);
        normZ /= length;
        normR /= length;

        if (i < NumElemFirstBound) {
            if (normZ < 0) {
                normZ = -normZ;
                normR = -normR;
            }
        } else if (i < NumElemFirstBound + NumElemSecondBound) {
            if (normR > 0) {
                normZ = -normZ;
                normR = -normR;
            }
        } else if (i < NumElemFirstBound + NumElemSecondBound + NumElemThirdBound) {
            if (normZ < 0) {
                normZ = -normZ;
                normR = -normR;
            }
        } else if (i < NumElemFirstBound + NumElemSecondBound + NumElemThirdBound + NumElemFourthBound) {
            if (normR < 0) {
                normZ = -normZ;
                normR = -normR;
            }
        } else if (i < NumElemFirstBound + NumElemSecondBound + NumElemThirdBound + NumElemFourthBound + NumElemFifthBound) {
            if (normZ > 0) {
                normZ = -normZ;
                normR = -normR;
            }
        } else if (i >= NumElemFirstBound + NumElemSecondBound + NumElemThirdBound + NumElemFourthBound + NumElemFifthBound) {
            if (normR > 0) {
                normZ = -normZ;
                normR = -normR;
            }
        }
        OutwardNormVec.push_back({normZ, normR});
    }
}

void BoundaryElementGeneration::buildSolidAngle() {
    SolidAngle.clear();
    SolidAngle.reserve(NumNodeTot);

    for (std::size_t i = 0; i < NumNodeTot; ++i) {
        Point2D TargetPoint = NodeCoord[i];
        Point2D PreviousPoint = NodeCoord[(i + NumNodeTot - 1) % NumNodeTot];
        Point2D NextPoint = NodeCoord[(i + 1) % NumNodeTot];
        
        double InternalAngle = calcInternalAngleOfTwoVertice(TargetPoint, PreviousPoint, NextPoint);
        SolidAngle.push_back(2 * InternalAngle);
    }
}

void BoundaryElementGeneration::dispNodeCoord() {
    std::cout << "BOUNDARY NODE COORIDNATES" << std::endl;
    for (std::size_t i = 0; i < NodeCoord.size(); ++i) {
        std::cout << "(" << NodeCoord[i].z << ", " << NodeCoord[i].r << "), " << std::endl;
    }
}

void BoundaryElementGeneration::dispOutwardNormVec() {
    std::cout << "OUTWARD NORMAL VECTOR" << std::endl;
    for (std::size_t i = 0; i < OutwardNormVec.size(); ++i) {
        std::cout << "(" << OutwardNormVec[i].nz << ", " << OutwardNormVec[i].nr << "), " << std::endl;
    }
}

void BoundaryElementGeneration::dispSolidAngle() {
    std::cout << "SOLID ANGLE OF EACH NODE" << std::endl;
    for (std::size_t i = 0; i < SolidAngle.size(); ++i) {
        std::cout << SolidAngle[i] << std::endl;
    }
}

void BoundaryElementGeneration::dispElemConn() {
    std::cout << "ELEMENT CONNECTION" << std::endl;
    for (std::size_t i = 0; i < ElemConn.size(); ++i) {
        std::cout << "[" << ElemConn[i].first << "---" << ElemConn[i].second << "]" << std::endl;
    }
}

void BoundaryElementGeneration::dispFlux() {
    std::cout << "FLUX" << std::endl;
    for (std::size_t i = 0; i < Flux.size(); ++i) {
        std::cout << "[" << Flux[i] << "]" << std::endl;
    }
}

void BoundaryElementGeneration::getNodeCoord(Point2DArr &NodeCoord) {
    NodeCoord = this->NodeCoord;
}

void BoundaryElementGeneration::getOutwardNormVec(NormVecArr &OutwardNormVec) {
    OutwardNormVec = this->OutwardNormVec;
}

void BoundaryElementGeneration::getSolidAngle(DoubleArr &SolidAngle) {
    SolidAngle = this->SolidAngle;
}

void BoundaryElementGeneration::getElemConn(ElemArr &ElemConn) {
    ElemConn = this->ElemConn;
}

void BoundaryElementGeneration::getFlux(DoubleArr &Flux) {
    Flux = this->Flux;
}

void BoundaryElementGeneration::initBoundaryElementGeneration(double L1, double L2, double H1, double H2, double H3, Size NumElemFirstBound, Size NumElemSecondBound, Size NumElemThirdBound, Size NumElemFourBound, Size NumElemFifthBound, Size NumElemSixthBound) {
    initPhysicalDomainDimension(L1, L2, H1, H2, H3);
    initBoundaryElementDimension(NumElemFirstBound, NumElemSecondBound, NumElemThirdBound, NumElemFourBound, NumElemFifthBound, NumElemSixthBound);
    buildElemConn();
    buildNodeCoord();
    buildOutwardNormVec();
    buildSolidAngle();
}

void BoundaryElementGeneration::initVelocity(double InVec) {
    double OutVec = InVec / (1.0 - std::pow(L2/L1, 2));
    Flux.clear();
    Flux.resize(NumNodeTot);
    for (std::size_t i = NumElemFirstBound + NumElemSecondBound + 1; i < NumElemFirstBound + NumElemSecondBound + NumElemThirdBound; ++i) {
        Flux[i] = OutVec;
    }
    for (std::size_t i = NumElemFirstBound + NumElemSecondBound + NumElemThirdBound + NumElemFourthBound + 1; i < NumElemFirstBound + NumElemSecondBound + NumElemThirdBound + NumElemFourthBound + NumElemFifthBound; ++i) {
        Flux[i] = -InVec;
    }
}
