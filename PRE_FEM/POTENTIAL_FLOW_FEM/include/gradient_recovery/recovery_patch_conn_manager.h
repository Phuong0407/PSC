#ifndef RECOVERY_PATCH_CONN_MANAGER_H
#define RECOVERY_PATCH_CONN_MANAGER_H

#include "../include/common_type_and_method.h"

#include <array>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <iostream>

class RecoveryPatchConnManager {
private:
    void genFirstLevelPatchNode(const Ind IndNode/*, const ElemConn &ElemConnData, NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnSetArr &FullNodeRecoveryConn*/);
    void genFirstLevelPatchNode(/*const Ind IndNode, const ElemConn &ElemConnData, NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnSetArr &FullNodeRecoveryConn*/);
    void findAllInternalNode(Ind IndNode, /*const NodeConnSetArr &FirstLevelNodeRecoveryConn, */NodeConnSet &InternalIndSet) const;
    // void updateNodeConnection(Ind IndNode, const NodeConnSet& InternalIndSet);
    void updateNodeConnection(const Ind IndNode, const NodeConnSet& InternalIndSet/*, const NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnSetArr &FullNodeRecoveryConn*/);
    void processFallbackConnection(Ind IndNode/*, const NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnSetArr &FullNodeRecoveryConn*/);

protected:
    const Size OrderOfPatch = 7;
    Size NumElemTol, NumNodeTol;
    ElemConn ElemConnData;
    NodeConnSetArr FullNodeRecoveryConn;
    NodeConnSetArr FirstLevelNodeRecoveryConn;
    NodeConnArr NodeRecoveryConnArr;

public:
    RecoveryPatchConnManager() = default;
    ~RecoveryPatchConnManager() { std::cout << "DESTROY RECOVERY PATCH CONNECTION MANAGER" << std::endl; }

public:
    // void initRecoveryPatchConnManager(Ind NumNodeTol, const ElemConn &ElemConn);
    void initRecoveryPatchConnManager(Size NumNodeTol, const ElemConn &ElemConnData/*, NodeConnSetArr &FullNodeRecoveryConn, NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnArr &NodeRecoveryConnArr*/);
    void genNodeRecoveryConn(/*Size NumNodeTol, const ElemConn &ElemConnData, NodeConnSetArr &FirstLevelNodeRecoveryConn, NodeConnSetArr &FullNodeRecoveryConn*/);
    void genNodeRecoveryConnArr(/*Size NumNodeTol, const NodeConnSetArr &FullNodeRecoveryConn, NodeConnArr &NodeRecoveryConnArr*/);
    void dispFullNodeRecoveryConn(/*const NodeConnSetArr &FullNodeRecoveryConn*/) const;
    void dispFirstLevelNodeRecoveryConn(/*const NodeConnSetArr &FirstLevelNodeRecoveryConn*/) const;
    void dispFullLevelNodeRecoveryConnArr(/*const NodeConnArr &NodeRecoveryConnArr*/) const;

public:
    // Size getSizeOfPatchOfCorrespondingNode(Ind IndNode) {
    //     return FullNodeRecoveryConn[IndNode].size();
    // }
    // void getFullNodeRecoveryConn(NodeConnArr &NodeRecoveryConnArr);
    // void getFirstLevelNodeRecoveryConn(NodeConnSetArr &FirstLevelNodeRecoveryConn);
};

#endif
