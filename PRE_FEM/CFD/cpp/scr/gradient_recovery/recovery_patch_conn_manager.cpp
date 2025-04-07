#include "../include/gradient_recovery/recovery_patch_conn_manager.h"

// void RecoveryPatchConnManager::initRecoveryPatchConnManager(Ind NumNodeTol, const ElemConn &ElemConnData) {
//     this->NumNodeTol = NumNodeTol;
//     this->NumElemTol = ElemConnData.size();
//     this->ElemConnData.assign(ElemConnData.begin(), ElemConnData.end());
//     FullNodeRecoveryConn.resize(NumNodeTol);
//     FirstLevelNodeRecoveryConn.resize(NumNodeTol);
//     genNodeRecoveryConn();
// }

void RecoveryPatchConnManager::initRecoveryPatchConnManager(Size NumNodeTol, const ElemConn &ElemConnData/*, NodeConnSetArr &FullNodeRecoveryConn, NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnArr &NodeRecoveryConnArr*/) {
    this->NumNodeTol = NumNodeTol;
    this->NumElemTol = ElemConnData.size();
    // this->ElemConnData = ElemConnData;
    this->ElemConnData.assign(ElemConnData.begin(), ElemConnData.end());
    this->FullNodeRecoveryConn.resize(NumNodeTol);
    this->FirstLevelNodeRecoveryConn.resize(NumNodeTol);
    this->NodeRecoveryConnArr.resize(NumNodeTol);
    // for (const auto& Elem : this->ElemConnData) {
        // std::cout << Elem[0] << " " << Elem[1] << " " << Elem[2] << std::endl;
    // }
    // std::cout << "END PRINT this->ElemConnData" << std::endl;

    genNodeRecoveryConn(/*NumNodeTol, ElemConnData, FirstLevelNodeRecoveryConn, FullNodeRecoveryConn*/);
    genNodeRecoveryConnArr(/*NumNodeTol, FullNodeRecoveryConn, NodeRecoveryConnArr*/);
    // dispFullLevelNodeRecoveryConnArr();
    // dispFirstLevelNodeRecoveryConn();
    // dispFullNodeRecoveryConn();
}

void RecoveryPatchConnManager::genNodeRecoveryConnArr(/*Size NumNodeTol, const NodeConnSetArr &FullNodeRecoveryConn, NodeConnArr &NodeRecoveryConnArr*/) {
    for (Ind NodeInd = 0; NodeInd < NumNodeTol; ++NodeInd)
        NodeRecoveryConnArr[NodeInd].push_back(NodeInd);

    for (Ind NodeInd = 0; NodeInd < NumNodeTol; ++NodeInd) {
        for (const auto &IndConnForNode : FullNodeRecoveryConn[NodeInd]) {
            if (IndConnForNode != NodeInd)
                NodeRecoveryConnArr[NodeInd].push_back(IndConnForNode);
            else
                continue;
        }
    }
}

void RecoveryPatchConnManager::genNodeRecoveryConn(/*Size NumNodeTol, const ElemConn &ElemConnData, NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnSetArr &FullNodeRecoveryConn*/) {
    // for (unsigned int IndNode = 0; IndNode < NumNodeTol; ++IndNode)
    //     genFirstLevelPatchNode(IndNode/*, ElemConnData, FirstLevelNodeRecoveryConn, FullNodeRecoveryConn*/);

    genFirstLevelPatchNode();
    // std::cout << "genFirstLevelPatchNode" << std::endl;

    for (unsigned int IndNode = 0; IndNode < NumNodeTol; ++IndNode) {
        if (FirstLevelNodeRecoveryConn[IndNode].size() < OrderOfPatch) {
            std::unordered_set<unsigned int> InternalIndSet;
            findAllInternalNode(IndNode, /*FirstLevelNodeRecoveryConn,*/ InternalIndSet);
            if (!InternalIndSet.empty()) {
                updateNodeConnection(IndNode, InternalIndSet/*, FirstLevelNodeRecoveryConn, FullNodeRecoveryConn*/);
            } else {
                processFallbackConnection(IndNode/*, FirstLevelNodeRecoveryConn, FullNodeRecoveryConn*/);
            }
        }
    }
}

void RecoveryPatchConnManager::genFirstLevelPatchNode(const Ind IndNode/*, const ElemConn &ElemConnData, NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnSetArr &FullNodeRecoveryConn*/) {
    FirstLevelNodeRecoveryConn[IndNode].insert({IndNode});
    for (const auto& Elem : ElemConnData) {
        if (std::find(Elem.begin(), Elem.end(), IndNode) != Elem.end()) {
            FirstLevelNodeRecoveryConn[IndNode].insert(Elem.begin(), Elem.end());
            FullNodeRecoveryConn[IndNode].insert(Elem.begin(), Elem.end());
        }
    }
}

void RecoveryPatchConnManager::genFirstLevelPatchNode(/*const Ind IndNode, const ElemConn &ElemConnData, NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnSetArr &FullNodeRecoveryConn*/) {
    for (Ind IndNode = 0; IndNode < NumNodeTol; ++IndNode) {
        FirstLevelNodeRecoveryConn[IndNode].insert({IndNode});
        for (const auto& Elem : ElemConnData) {
            if (std::find(Elem.begin(), Elem.end(), IndNode) != Elem.end()) {
                FirstLevelNodeRecoveryConn[IndNode].insert(Elem.begin(), Elem.end());
                FullNodeRecoveryConn[IndNode].insert(Elem.begin(), Elem.end());
            }
        }
    }
}

void RecoveryPatchConnManager::findAllInternalNode(const Ind IndNode, /*const NodeConnSetArr &FirstLevelNodeRecoveryConn, */NodeConnSet &InternalIndSet) const {
    for (std::size_t IndConnectedNode : FirstLevelNodeRecoveryConn[IndNode]) {
        if (IndConnectedNode != IndNode) {
            if (FirstLevelNodeRecoveryConn[IndConnectedNode].size() >= OrderOfPatch) {
                InternalIndSet.insert(IndConnectedNode);
            }
        }
    }
}

void RecoveryPatchConnManager::updateNodeConnection(const Ind IndNode, const NodeConnSet &InternalIndSet/*, const NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnSetArr &FullNodeRecoveryConn*/) {
    for (const auto& InternalInd : InternalIndSet) {
        FullNodeRecoveryConn[IndNode].insert(FirstLevelNodeRecoveryConn[InternalInd].begin(), FirstLevelNodeRecoveryConn[InternalInd].end());
    }
}

void RecoveryPatchConnManager::processFallbackConnection(Ind IndNode) {
    for (const auto& IndConnectedNode : FirstLevelNodeRecoveryConn[IndNode]) {
        if (IndConnectedNode != IndNode) {
            std::unordered_set<unsigned int> TempInternalIndSet;
            findAllInternalNode(IndConnectedNode, /*FirstLevelNodeRecoveryConn, */TempInternalIndSet);
            
            if (!TempInternalIndSet.empty()) {
                updateNodeConnection(IndNode, TempInternalIndSet/*, FirstLevelNodeRecoveryConn, FullNodeRecoveryConn*/);
                if (FullNodeRecoveryConn.size() >= OrderOfPatch) {
                    return;
                }
            }
        }
    }
}

void RecoveryPatchConnManager::dispFirstLevelNodeRecoveryConn(/*const NodeConnSetArr &FirstLevelNodeRecoveryConn*/) const {
    std::cout << "FIRST LEVEL NODE RECOVERY CONNECTION" << std::endl;
    for (size_t i = 0; i < FirstLevelNodeRecoveryConn.size(); ++i) {
        std::cout << "Node " << i << ": [ ";
        for (const auto& ConnectedNode : FirstLevelNodeRecoveryConn[i]) {
            std::cout << ConnectedNode << " ";
        }
        std::cout << "]" << std::endl;
    }
}

void RecoveryPatchConnManager::dispFullNodeRecoveryConn(/*const NodeConnSetArr &FullNodeRecoveryConn*/) const {
    std::cout << "FULL NODE RECOVERY CONNECTION" << std::endl;
    for (size_t i = 0; i < FullNodeRecoveryConn.size(); ++i) {
        std::cout << "Node " << i << ": [ ";
        for (const auto& ConnectedNode : FullNodeRecoveryConn[i]) {
            std::cout << ConnectedNode << " ";
        }
        std::cout << "]" << std::endl;
    }
}

void RecoveryPatchConnManager::dispFullLevelNodeRecoveryConnArr(/*const NodeConnArr &NodeRecoveryConnArr*/) const {
    std::cout << "FULL NODE RECOVERY CONNECTION ARRAY" << std::endl;
    for (size_t i = 0; i < NodeRecoveryConnArr.size(); ++i) {
        std::cout << "Node " << i << ": [ ";
        for (const auto& ConnectedNode : NodeRecoveryConnArr[i]) {
            std::cout << ConnectedNode << " ";
        }
        std::cout << "]" << std::endl;
    }
}

// void RecoveryPatchConnManager::getFullNodeRecoveryConn(NodeConnArr &FullNodeRecoveryConnArr) {
//     FullNodeRecoveryConnArr.assign(this->NodeRecoveryConnArr.begin(), this->NodeRecoveryConnArr.end());
// }

// void RecoveryPatchConnManager::getFirstLevelNodeRecoveryConn(NodeConnSetArr &FirstLevelNodeRecoveryConn) {
//     FirstLevelNodeRecoveryConn.assign(this->FirstLevelNodeRecoveryConn.begin(), this->FirstLevelNodeRecoveryConn.end());
// }
