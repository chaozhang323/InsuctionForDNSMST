#include "PartitionTecUtilDecorator.h"
namespace tecplot { namespace ___3933 { PartitionTecUtilDecorator::PartitionTecUtilDecorator( ___37& ___36, ___4636 zone) : ___2337(___36) , m_zoneNum(zone) { REQUIRE(___36.___896()); REQUIRE(0 < zone && ___36.___4638(zone)); REQUIRE(___36.zoneIsPartitioned(zone)); } PartitionTecUtilDecorator::~PartitionTecUtilDecorator() {} void PartitionTecUtilDecorator::___3817(char** ___3855) const { ___2337.___3817(___3855); } void PartitionTecUtilDecorator::___3827(___3839* ___3819) const { ___2337.___3827(___3819); } ___2227 PartitionTecUtilDecorator::___3832(___3839 ___3819) const { return ___2337.___3832(___3819); } char* PartitionTecUtilDecorator::___3833(___3839 ___3819, ___2227 ___3853) const { return ___2337.___3833(___3819, ___3853); } ___264 PartitionTecUtilDecorator::___235() const { return ___2337.___235(); } ___264 PartitionTecUtilDecorator::___274(___4636  ) const { ___478(___1305); return ___2337.___274(m_zoneNum); } ___264 PartitionTecUtilDecorator::___273(___4352 ___4336) const { return ___2337.___273(___4336); } int32_t PartitionTecUtilDecorator::___247(___264 ___265) const { return ___2337.___247(___265); } void PartitionTecUtilDecorator::___243(___264 ___265, int32_t index, char** ___2685, ___90* ___4314, AuxDataType_e* type, ___372* ___3361) const { ___2337.___243(___265, index, ___2685, ___4314, type, ___3361); } ___372 PartitionTecUtilDecorator::___896(void) const { return ___2337.___896(); } ___372 PartitionTecUtilDecorator::datasetGetTitle(char** datasetTitle) const { return ___2337.datasetGetTitle(datasetTitle); } int32_t PartitionTecUtilDecorator::___889(void) const { return ___2337.___889(); } ___3501 PartitionTecUtilDecorator::datasetGetRelevantZones(double solutionTimeMin, double solutionTimeMax, ___372 ignoreStaticZones) const { return ___2337.datasetGetRelevantZones(solutionTimeMin, solutionTimeMax, ignoreStaticZones); } ___4636 PartitionTecUtilDecorator::___891(void) const { return ___2337.zoneGetNumPartitions(m_zoneNum); } ___4352 PartitionTecUtilDecorator::___890(void) const { return ___2337.___890(); } ___4352 PartitionTecUtilDecorator::___4345(char ___214) const { return ___2337.___4345(___214); } ___372 PartitionTecUtilDecorator::___4344(___4352 ___4368, char** ___4362) const { return ___2337.___4344(___4368, ___4362); } int32_t PartitionTecUtilDecorator::___4343(___4352 ___4368) const { return ___2337.___4343(___4368); } ___372 PartitionTecUtilDecorator::___4638(___4636  ) const { return ___2337.___4638(m_zoneNum); } ___372 PartitionTecUtilDecorator::___4640(___4636  ) const { return ___2337.___4640(m_zoneNum); } ___372 PartitionTecUtilDecorator::___4641(___4636 ___4658) const { ___4278(___4658); return ___2337.___4641(m_zoneNum); } int32_t PartitionTecUtilDecorator::___4613(___4636  ) const { return ___2337.___4613(m_zoneNum); } ___372 PartitionTecUtilDecorator::___4614(___3501* ___1153) const { return ___2337.___4614(___1153); } void PartitionTecUtilDecorator::___4615(___4636 partitionNum, ___1844& ___2715) const { ___2337.zonePartitionGetIJK(m_zoneNum, partitionNum, ___2715); } ___372 PartitionTecUtilDecorator::___4616(___4636  , char** ___4652) const { return ___2337.___4616(m_zoneNum, ___4652); } ___4636 PartitionTecUtilDecorator::___4617(___4636  ) const { return ___2337.___4617(m_zoneNum); } double PartitionTecUtilDecorator::___4618(___4636  ) const { return ___2337.___4618(m_zoneNum); } ___1172 PartitionTecUtilDecorator::___4619(___4636  ) const { return ___2337.___4619(m_zoneNum); } ZoneType_e PartitionTecUtilDecorator::___4620(___4636  ) const { return ___2337.___4620(m_zoneNum); } ___372 PartitionTecUtilDecorator::___4353(___4352 ___4368) const { return ___2337.___4353(___4368); } ___372 PartitionTecUtilDecorator::varGetEnabled(___3501* enabledVars) const { return ___2337.varGetEnabled(enabledVars); } int32_t PartitionTecUtilDecorator::solutionTimeGetNumTimeSteps() const { return ___2337.solutionTimeGetNumTimeSteps(); } double PartitionTecUtilDecorator::solutionTimeGetMinTime() const { return ___2337.solutionTimeGetMinTime();
} double PartitionTecUtilDecorator::solutionTimeGetMaxTime() const { return ___2337.solutionTimeGetMaxTime(); } ___372 PartitionTecUtilDecorator::___3768() const { return ___2337.___3768(); } GeomID_t PartitionTecUtilDecorator::___1592(void) { return ___2337.___1592(); } TextID_t PartitionTecUtilDecorator::___4087(void) { return ___2337.___4087(); } int32_t PartitionTecUtilDecorator::___797(void) { return ___2337.___797(); } ___372 PartitionTecUtilDecorator::___796(___3839* ___2171, int32_t ___4453) { return ___2337.___796(___2171, ___4453); } void PartitionTecUtilDecorator::___3779(char const* ___3001, ___372 ___3584, ___372 ___3579) const { ___2337.___3779(___3001, ___3584, ___3579); } void PartitionTecUtilDecorator::___3778(char const* ___3001) const { ___2337.___3778(___3001); } ___372 PartitionTecUtilDecorator::___3769(int PercentDone) const { return ___2337.___3769(PercentDone); } void PartitionTecUtilDecorator::___3770(void) const { ___2337.___3770(); } ___372 PartitionTecUtilDecorator::___1983(void) const { return ___2337.___1983(); } void PartitionTecUtilDecorator::___858(void) { ___2337.___858(); } void PartitionTecUtilDecorator::___859(void) { ___2337.___859(); } ___4636 PartitionTecUtilDecorator::___544(___3501 ___4684, ___4636  ) const { return ___2337.___544(___4684, m_zoneNum); } ___3501 PartitionTecUtilDecorator::___545(___4636  ) const { return ___2337.___545(m_zoneNum); } ValueLocation_e PartitionTecUtilDecorator::___910(___4636  , ___4352 ___4368) const { return ___2337.___910(m_zoneNum, ___4368); } ValueLocation_e PartitionTecUtilDecorator::___911(___1361 ___1351) const { return ___2337.___911(___1351); } ___372 PartitionTecUtilDecorator::___913(___4636 partitionNum, ___4352 ___4336, double* minVal, double* maxVal) const { return ___2337.dataValueGetMinMaxByZonePartitionVar(m_zoneNum, partitionNum, ___4336, minVal, maxVal); } ___372 PartitionTecUtilDecorator::___912(___1361 ___1351, double* minVal, double* maxVal) const { return ___2337.___912(___1351, minVal, maxVal); } FieldDataType_e PartitionTecUtilDecorator::___923(___4636  , ___4352 ___4368) { return ___2337.___923(m_zoneNum, ___4368); } ___1172 PartitionTecUtilDecorator::___921(___4636  , ___4352 ___4368) { return ___2337.___921(m_zoneNum, ___4368); } ___1361 PartitionTecUtilDecorator::___918(___4636 partitionNum, ___4352 ___4368) { return ___2337.dataValuePartitionGetReadableNLRef(m_zoneNum, partitionNum, ___4368); } ___1361 PartitionTecUtilDecorator::___915(___4636 partitionNum, ___4352 ___4368) { return ___2337.dataValuePartitionGetReadableCCRef(m_zoneNum, partitionNum, ___4368); } ___1361 PartitionTecUtilDecorator::___917(___4636 partitionNum, ___4352 ___4368) { return ___2337.dataValuePartitionGetReadableNativeRef(m_zoneNum, partitionNum, ___4368); } ___1361 PartitionTecUtilDecorator::___916(___4636 partitionNum, ___4352 ___4368) { return ___2337.dataValuePartitionGetReadableDerivedRef(m_zoneNum, partitionNum, ___4368); } ___1361 PartitionTecUtilDecorator::___924(___4636 partitionNum, ___4352 ___4368) { return ___2337.dataValuePartitionGetWritableNativeRef(m_zoneNum, partitionNum, ___4368); } double PartitionTecUtilDecorator::___909(___1361 ___1351, ___81 ___2733) { return ___2337.___909(___1351, ___2733); } void PartitionTecUtilDecorator::dataValueSetByRef(___1361 ___1351, ___81 ___2733, double ___4298) { ___2337.dataValueSetByRef(___1351, ___2733, ___4298); } void PartitionTecUtilDecorator::___919(___4636 partitionNum, ___4352 ___4368, void** ___880, FieldDataType_e* ___1363) { ___2337.dataValuePartitionGetReadableRawPtr(m_zoneNum, partitionNum, ___4368, ___880, ___1363); } void PartitionTecUtilDecorator::___925(___4636 partitionNum, ___4352 ___4368, void** ___880, FieldDataType_e* ___1363) { ___2337.dataValuePartitionGetWritableRawPtr(m_zoneNum, partitionNum, ___4368, ___880, ___1363); } ___4636 PartitionTecUtilDecorator::___914(___3501 ___4684, ___4636  , ___4352 ___4336) const
{ return ___2337.___914(___4684, m_zoneNum, ___4336); } ___3501 PartitionTecUtilDecorator::___922(___4636  , ___4352 ___4368) const { return ___2337.___922(m_zoneNum, ___4368); } ___372 PartitionTecUtilDecorator::___926(___4636  , ___4352 ___4368) const { return ___2337.___926(m_zoneNum, ___4368); } ___1383 PartitionTecUtilDecorator::___927(___1361 ___1309) { return ___2337.___927(___1309); } ___1384 PartitionTecUtilDecorator::___928(___1361 ___1309) { return ___2337.___928(___1309); } FieldDataType_e PartitionTecUtilDecorator::___920(___1361 ___1352) { return ___2337.___920(___1352); } ___2727 PartitionTecUtilDecorator::___867(___4636 partitionNum) { return ___2337.dataNodePartitionGetReadableRef(m_zoneNum, partitionNum); } ___2727 PartitionTecUtilDecorator::___869(___4636 partitionNum) { return ___2337.dataNodePartitionGetWritableRef(m_zoneNum, partitionNum); } ___2718 PartitionTecUtilDecorator::___865(___2727 ___2723, ___465 ___468, ___682 ___683) { return ___2337.___865(___2723, ___468, ___683); } void PartitionTecUtilDecorator::___870(___2727 ___2723, ___465 ___468, ___682 ___683, ___2718 ___2716) { ___2337.___870(___2723, ___468, ___683, ___2716); } OffsetDataType_e PartitionTecUtilDecorator::dataNodeGetRawItemType(___2727 ___2723) { return ___2337.dataNodeGetRawItemType(___2723); } int32_t* PartitionTecUtilDecorator::dataNodeGetRawPtrByRef(___2727 ___2723) { return ___2337.dataNodeGetRawPtrByRef(___2723); } int64_t* PartitionTecUtilDecorator::dataNodeGetRawPtrByRef64(___2727 ___2723) { return ___2337.dataNodeGetRawPtrByRef64(___2723); } ___2742 PartitionTecUtilDecorator::dataNodeToElemMapGetReadableRef(___4636 partitionNum) const { return ___2337.dataNodeToElemMapPartitionGetReadableRef(m_zoneNum, partitionNum); } ___465 PartitionTecUtilDecorator::dataNodeToElemMapGetNumElems(___2742 nodeToElemMap, ___2718 ___2709) const { return ___2337.dataNodeToElemMapGetNumElems(nodeToElemMap, ___2709); } ___465 PartitionTecUtilDecorator::dataNodeToElemMapGetElem(___2742 nodeToElemMap, ___2718 ___2709, ___465 elemOffset) const { return ___2337.dataNodeToElemMapGetElem(nodeToElemMap, ___2709, elemOffset); } FaceNeighborMode_e PartitionTecUtilDecorator::___836(___4636  ) const { return ___2337.___836(m_zoneNum); } void PartitionTecUtilDecorator::___837(___1292 ___1274, ___2227 ___1144, int32_t face, int32_t ___2692, ___2227* ___2691, ___4636* ___2695) const { ___2337.___837(___1274, ___1144, face, ___2692, ___2691, ___2695); } ___372 PartitionTecUtilDecorator::___835(___1292 ___1274, ___2227 ___1144, int32_t face, ___3501 ___4) const { return ___2337.___835(___1274, ___1144, face, ___4); } int32_t PartitionTecUtilDecorator::___838(___1292 ___1274, ___2227 ___1144, int32_t face, ___372* neighborsAreUserSpecified) const { return ___2337.___838(___1274, ___1144, face, neighborsAreUserSpecified); } ___1292 PartitionTecUtilDecorator::___839(___4636  ) const { return ___2337.___839(m_zoneNum); } void PartitionTecUtilDecorator::___3484(___3501* set) const { ___2337.___3484(set); } ___3493 PartitionTecUtilDecorator::___3491(___3501 set, ___3493 ___2401) const { return ___2337.___3491(set, ___2401); } ___3493 PartitionTecUtilDecorator::setGetPrevMember(___3501 set, ___3493 ___2401) const { return ___2337.setGetPrevMember(set, ___2401); } ___3493 PartitionTecUtilDecorator::setGetMemberCount(___3501 set) const { return ___2337.setGetMemberCount(set); } ___372 PartitionTecUtilDecorator::___3495(___3501 set, ___3493 ___2401) const { return ___2337.___3495(set, ___2401); } ___372 PartitionTecUtilDecorator::setIsEqual(___3501 ___3477, ___3501 ___3478) const { return ___2337.setIsEqual(___3477, ___3478); } void PartitionTecUtilDecorator::setRemoveMember(___3501 set, ___3493 ___2401) const { ___2337.setRemoveMember(set, ___2401); } void PartitionTecUtilDecorator::___1557(GeomID_t ___1805, int32_t ___3157, ___2227 ___3141, double* x, double* ___4583) const { ___2337.___1557(___1805, ___3157, ___3141, x, ___4583); } void PartitionTecUtilDecorator::___1558(GeomID_t ___1805, ___2227 ___3141, double* x, double* ___4583) const { ___2337.___1558(___1805, ___3141, x, ___4583);
} void PartitionTecUtilDecorator::___1560(GeomID_t ___1805, int32_t ___3157, ___2227 ___3141, double* x, double* ___4583, double* z) const { ___2337.___1560(___1805, ___3157, ___3141, x, ___4583, z); } void PartitionTecUtilDecorator::___1561(GeomID_t ___1805, ___2227 ___3141, double* x, double* ___4583, double* z) const { ___2337.___1561(___1805, ___3141, x, ___4583, z); } double PartitionTecUtilDecorator::___1564(GeomID_t ___1805) const { return ___2337.___1564(___1805); } ArrowheadAttachment_e PartitionTecUtilDecorator::___1565(GeomID_t ___1805) const { return ___2337.___1565(___1805); } double PartitionTecUtilDecorator::___1566(GeomID_t ___1805) const { return ___2337.___1566(___1805); } ArrowheadStyle_e PartitionTecUtilDecorator::___1567(GeomID_t ___1805) const { return ___2337.___1567(___1805); } double PartitionTecUtilDecorator::___1570(GeomID_t ___1805) const { return ___2337.___1570(___1805); } int32_t PartitionTecUtilDecorator::___1576(GeomID_t ___1805) const { return ___2337.___1576(___1805); } void PartitionTecUtilDecorator::___1577(GeomID_t ___1805, double* ___1824, double* ___4394) const { ___2337.___1577(___1805, ___1824, ___4394); } void PartitionTecUtilDecorator::___1591(GeomID_t ___1805, double* ___4574, double* ___4591, double* ___4715) const { return ___2337.___1591(___1805, ___4574, ___4591, ___4715); } Clipping_e PartitionTecUtilDecorator::___1593(GeomID_t ___1805) const { return ___2337.___1593(___1805); } ___516 PartitionTecUtilDecorator::___1594(GeomID_t ___1805) const { return ___2337.___1594(___1805); } DrawOrder_e PartitionTecUtilDecorator::___1595(GeomID_t ___1805) const { return ___2337.___1595(___1805); } ___516 PartitionTecUtilDecorator::___1596(GeomID_t ___1805) const { return ___2337.___1596(___1805); } ___372 PartitionTecUtilDecorator::___1597(GeomID_t ___1805) const { return ___2337.___1597(___1805); } LinePattern_e PartitionTecUtilDecorator::___1598(GeomID_t ___1805) const { return ___2337.___1598(___1805); } double PartitionTecUtilDecorator::___1599(GeomID_t ___1805) const { return ___2337.___1599(___1805); } ___372 PartitionTecUtilDecorator::___1600(GeomID_t ___1805, char** macroFunctionCmd) const { return ___2337.___1600(___1805, macroFunctionCmd); } GeomID_t PartitionTecUtilDecorator::___1601(GeomID_t ___1805) const { return ___2337.___1601(___1805); } double PartitionTecUtilDecorator::___1602(GeomID_t ___1805) const { return ___2337.___1602(___1805); } CoordSys_e PartitionTecUtilDecorator::___1603(GeomID_t ___1805) const { return ___2337.___1603(___1805); } GeomID_t PartitionTecUtilDecorator::___1604(GeomID_t ___1805) const { return ___2337.___1604(___1805); } Scope_e PartitionTecUtilDecorator::___1605(GeomID_t ___1805) const { return ___2337.___1605(___1805); } GeomForm_e PartitionTecUtilDecorator::___1606(GeomID_t ___1805) const { return ___2337.___1606(___1805); } ___4636 PartitionTecUtilDecorator::___1607(GeomID_t ___1805) const { return ___2337.___1607(___1805); } ___372 PartitionTecUtilDecorator::___1610(GeomID_t ___1805) const { return ___2337.___1610(___1805); } ___2227 PartitionTecUtilDecorator::___1619(GeomID_t ___1805, int32_t ___3157) const { return ___2337.___1619(___1805, ___3157); } ___2227 PartitionTecUtilDecorator::___1620(GeomID_t ___1805) const { return ___2337.___1620(___1805); } ___2227 PartitionTecUtilDecorator::___1626(GeomID_t ___1805) const { return ___2337.___1626(___1805); } void PartitionTecUtilDecorator::___1628(GeomID_t ___1805, double* ___4458, double* ___1826) const { return ___2337.___1628(___1805, ___4458, ___1826); } double PartitionTecUtilDecorator::___1648(GeomID_t ___1805) const { return ___2337.___1648(___1805); } ___516 PartitionTecUtilDecorator::___4064(TextID_t ___4171) const { return ___2337.___4064(___4171); } ___516 PartitionTecUtilDecorator::___4065(TextID_t ___4171) const { return ___2337.___4065(___4171); } double PartitionTecUtilDecorator::___4066(TextID_t ___4171) const { return ___2337.___4066(___4171); } double PartitionTecUtilDecorator::___4067(TextID_t ___4171) const { return ___2337.___4067(___4171); } TextBox_e PartitionTecUtilDecorator::___4068(TextID_t ___4171) const { return ___2337.___4068(___4171); } TextAnchor_e PartitionTecUtilDecorator::___4084(TextID_t ___4171) const { return ___2337.___4084(___4171); } void PartitionTecUtilDecorator::___4085(TextID_t ___4171, double* ___4574, double* ___4591, double* ___4715) const
{ return ___2337.___4085(___4171, ___4574, ___4591, ___4715); } double PartitionTecUtilDecorator::___4086(TextID_t ___4171) const { return ___2337.___4086(___4171); } Clipping_e PartitionTecUtilDecorator::___4088(TextID_t ___4171) const { return ___2337.___4088(___4171); } ___516 PartitionTecUtilDecorator::___4089(TextID_t ___4171) const { return ___2337.___4089(___4171); } double PartitionTecUtilDecorator::___4090(TextID_t ___4171) const { return ___2337.___4090(___4171); } double PartitionTecUtilDecorator::___4091(TextID_t ___4171) const { return ___2337.___4091(___4171); } ___372 PartitionTecUtilDecorator::___4092(TextID_t ___4171, char** ___2330) const { return ___2337.___4092(___4171, ___2330); } TextID_t PartitionTecUtilDecorator::___4093(TextID_t ___4171) const { return ___2337.___4093(___4171); } CoordSys_e PartitionTecUtilDecorator::___4094(TextID_t ___4171) const { return ___2337.___4094(___4171); } TextID_t PartitionTecUtilDecorator::___4095(TextID_t ___4171) const { return ___2337.___4095(___4171); } Scope_e PartitionTecUtilDecorator::___4096(TextID_t ___4171) const { return ___2337.___4096(___4171); } Units_e PartitionTecUtilDecorator::___4097(TextID_t ___4171) const { return ___2337.___4097(___4171); } ___372 PartitionTecUtilDecorator::___4098(TextID_t ___4171, char** ___4126) const { return ___2337.___4098(___4171, ___4126); } char* PartitionTecUtilDecorator::___4099(TextID_t ___4171) const { return ___2337.___4099(___4171); } ___372 PartitionTecUtilDecorator::___4100(TextID_t ___4171) const { return ___2337.___4100(___4171); } ___372 PartitionTecUtilDecorator::___4101(TextID_t ___4171) const { return ___2337.___4101(___4171); } ___4636 PartitionTecUtilDecorator::___4102(TextID_t ___4171) const { return ___2337.___4102(___4171); } ___372 PartitionTecUtilDecorator::___4105(TextID_t ___4171) const { return ___2337.___4105(___4171); } ___2664 PartitionTecUtilDecorator::___4152() { return ___2337.___4152(); } void PartitionTecUtilDecorator::___4153(___2664* mutex) { ___2337.___4153(mutex); } void PartitionTecUtilDecorator::___4154(___2664 mutex) { ___2337.___4154(mutex); } void PartitionTecUtilDecorator::___4155(___2664 mutex) { ___2337.___4155(mutex); } void PartitionTecUtilDecorator::___4156(___4160 ___2118, ___90 ___2123, ___2120 ___2119) { ___2337.___4156(___2118, ___2123, ___2119); } int PartitionTecUtilDecorator::___4157() { return ___2337.___4157(); } ___2120 PartitionTecUtilDecorator::___4158() { return ___2337.___4158(); } void PartitionTecUtilDecorator::___4159(___2120* ___2119) { ___2337.___4159(___2119); } void PartitionTecUtilDecorator::___4161(___2120 ___2119) { ___2337.___4161(___2119); } ___372 PartitionTecUtilDecorator::___4304(___1361 ___1351) const { return ___2337.___4304(___1351); } ___372 PartitionTecUtilDecorator::___4309(___2727 ___2723) const { return ___2337.___4309(___2723); } PlotType_e PartitionTecUtilDecorator::___1513() const { return ___2337.___1513(); } int32_t PartitionTecUtilDecorator::datasetGetNumPartitionFiles() const { return ___2337.datasetGetNumPartitionFiles(); } int32_t PartitionTecUtilDecorator::zoneGetOwnerProcess(___4636 partitionNum) const { return ___2337.zonePartitionGetOwnerProcess(m_zoneNum, partitionNum); } int32_t PartitionTecUtilDecorator::zonePartitionGetOwnerProcess(___4636  , ___4636 /*partitionNum*/) const { ___478(___1305); return 0; } ___372 PartitionTecUtilDecorator::zoneIsPartitioned(___4636  ) const { return ___1305; } ___4636 PartitionTecUtilDecorator::zoneGetNumPartitions(___4636  ) const { return 0; } void PartitionTecUtilDecorator::zonePartitionGetIJK(___4636  , ___4636 /*partitionNum*/, ___1844& /*___1861*/) const { ___478(___1305); } void PartitionTecUtilDecorator::zonePartitionGetIJKOffset(___4636 ___4658, ___4636 partitionNum, ___1844& ___1862) const { ___4278(___4658); ___4278(partitionNum); ___4278(___1862); ___478(___1305); } ___372 PartitionTecUtilDecorator::dataValueGetMinMaxByZonePartitionVar( ___4636  , ___4636  , ___4352  , double*  , double*  ) const { ___478(___1305); return ___1305; } ___1361 PartitionTecUtilDecorator::dataValuePartitionGetReadableNLRef( ___4636  , ___4636  , ___4352  ) const { ___478(___1305); return NULL; } ___1361 PartitionTecUtilDecorator::dataValuePartitionGetReadableCCRef(
___4636  , ___4636  , ___4352  ) const { ___478(___1305); return NULL; } ___1361 PartitionTecUtilDecorator::dataValuePartitionGetReadableNativeRef( ___4636  , ___4636  , ___4352  ) const { ___478(___1305); return NULL; } ___1361 PartitionTecUtilDecorator::dataValuePartitionGetReadableDerivedRef( ___4636  , ___4636  , ___4352  ) const { ___478(___1305); return NULL; } ___1361 PartitionTecUtilDecorator::dataValuePartitionGetWritableNativeRef( ___4636  , ___4636  , ___4352  ) const { ___478(___1305); return NULL; } void PartitionTecUtilDecorator::dataValuePartitionGetReadableRawPtr( ___4636  , ___4636  , ___4352  , void** ___880, FieldDataType_e* ___1363) const { ___478(___1305); *___880 = NULL; *___1363 = ___1369; } void PartitionTecUtilDecorator::dataValuePartitionGetWritableRawPtr( ___4636  , ___4636  , ___4352  , void** ___880, FieldDataType_e* ___1363) const { ___478(___1305); *___880 = NULL; *___1363 = ___1369; } ___2727 PartitionTecUtilDecorator::dataNodePartitionGetReadableRef( ___4636  , ___4636  ) const { ___478(___1305); return NULL; } ___2727 PartitionTecUtilDecorator::dataNodePartitionGetWritableRef( ___4636  , ___4636  ) const { ___478(___1305); return NULL; } ___2742 PartitionTecUtilDecorator::dataNodeToElemMapPartitionGetReadableRef( ___4636  , ___4636  ) const { ___478(___1305); return NULL; } GhostInfo_pa PartitionTecUtilDecorator::zoneGhostNodeInfoGetRef(___4636 partitionNum) const { return ___2337.zonePartitionGhostNodeInfoGetRef(m_zoneNum, partitionNum); } GhostInfo_pa PartitionTecUtilDecorator::zoneGhostCellInfoGetRef(___4636 partitionNum) const { return ___2337.zonePartitionGhostCellInfoGetRef(m_zoneNum, partitionNum); } GhostInfo_pa PartitionTecUtilDecorator::zonePartitionGhostNodeInfoGetRef( ___4636  , ___4636  ) const { ___478(___1305); return 0; } GhostInfo_pa PartitionTecUtilDecorator::zonePartitionGhostCellInfoGetRef( ___4636  , ___4636  ) const { ___478(___1305); return 0; } ___81 PartitionTecUtilDecorator::ghostInfoGetNumItemsByRef(GhostInfo_pa ghostInfo) const { return ___2337.ghostInfoGetNumItemsByRef(ghostInfo); } ___81 PartitionTecUtilDecorator::ghostInfoGetItemByRef(GhostInfo_pa ghostInfo, ___81 itemNum) const { return ___2337.ghostInfoGetItemByRef(ghostInfo, itemNum); } ___2090::___2980 PartitionTecUtilDecorator::ghostInfoGetNeighborByRef(GhostInfo_pa ghostInfo, ___81 itemNum) const { return ___2337.ghostInfoGetNeighborByRef(ghostInfo, itemNum); } ___81 PartitionTecUtilDecorator::ghostInfoGetNeighborItemByRef(GhostInfo_pa ghostInfo, ___81 itemNum) const { return ___2337.ghostInfoGetNeighborItemByRef(ghostInfo, itemNum); } }}
