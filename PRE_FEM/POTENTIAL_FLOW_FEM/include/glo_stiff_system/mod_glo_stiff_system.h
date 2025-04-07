#ifndef MODIFICATION_GLOBAL_STIFF_SYSTEM_H
#define MODIFICATION_GLOBAL_STIFF_SYSTEM_H

#include "../include/common_type_and_method.h"
#include <cmath>
#include <vector>
#include <iostream>

class ModificationGloStiffSystem {
private:
    inline void eliminateRowDirichletBound(const Ind RowInd, Mat &StiffMat) {
        std::fill(StiffMat[RowInd].begin(), StiffMat[RowInd].end(), 0.0);
    }
    inline void eliminateColDirichletBound(const Ind ColInd, Mat &StiffMat) {
        for (std::size_t j = 0; j < StiffMat.size(); ++j) {
            StiffMat[j][ColInd] = 0.0;
        }
    }
    inline void modifyRowColGloStiffMat(const Ind BoundNodeInd, Mat &StiffMat) {
        eliminateColDirichletBound(BoundNodeInd, StiffMat);
        eliminateRowDirichletBound(BoundNodeInd, StiffMat);
        StiffMat[BoundNodeInd][BoundNodeInd] = 1.0;
    }

    void modifyCorrespodingColFromDirichletBoundVal(const Ind IndCorrespondingCol, const double ValDirichletBound, const Mat &StiffMat, ColVect &ForceMat);

protected:
    ModificationGloStiffSystem() = default;
    void modifyGloStiffSystem(IndArr &IndDirichletBound, const ColVect &ValDirichletBound, Mat &StiffMat, ColVect &ForceMat);
    void dispModifiedForceMat(ColVect &ForceMat) const;
};

#endif