#include "../include/glo_stiff_system/mod_glo_stiff_system.h"

void ModificationGloStiffSystem::modifyGloStiffSystem(IndArr &IndDirichletBound, const ColVect &ValDirichletBound, Mat &StiffMat, ColVect &ForceMat) {
    for (Ind i = 0; i < IndDirichletBound.size(); ++i) {
        Ind CorrespondingColInd = IndDirichletBound[i];
        modifyCorrespodingColFromDirichletBoundVal(CorrespondingColInd, ValDirichletBound[i], StiffMat, ForceMat);
        modifyRowColGloStiffMat(CorrespondingColInd, StiffMat);
        ForceMat[CorrespondingColInd] = ValDirichletBound[i];
    }
}

void ModificationGloStiffSystem::dispModifiedForceMat(ColVect &ForceMat) const {
    std::cout << "MODIFIED FORCE MATRIX" << std::endl;
    for(std::size_t i = 0; i < ForceMat.size(); ++i) {
        std::cout << ForceMat[i] << std::endl;
    }
}

void ModificationGloStiffSystem::modifyCorrespodingColFromDirichletBoundVal(const Ind CorrespondingColInd, const double ValDirichletBound, const Mat &StiffMat, ColVect &ForceMat) {
    for (std::size_t j = 0; j < ForceMat.size(); ++j) {
        ForceMat[j] -= StiffMat[j][CorrespondingColInd] * ValDirichletBound;
    }
}
