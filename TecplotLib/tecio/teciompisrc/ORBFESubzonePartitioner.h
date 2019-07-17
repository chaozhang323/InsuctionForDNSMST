 #pragma once
#include "FESubzonePartitionerInterface.h"
#include "OrthogonalBisection.h"
namespace tecplot { namespace ___3933 { class OrbFESubzonePartitioner : public virtual FESubzonePartitionerInterface { UNCOPYABLE_CLASS(OrbFESubzonePartitioner); private: ___4636 const               ___2677; ___2090::ItemOffset_t const m_fixedSubzoneSize; OrthogonalBisection             m_cellOrb; OrthogonalBisection             m_nodeOrb; ItemAddressArray                m_szCoordsOfOrginalZoneCells; ItemAddressArray                m_szCoordsOfOrginalZoneNodes; ___372 partitionIntoSubzones(___37& ___36); public: OrbFESubzonePartitioner( ___37&               ___36, ___4636               zone, ___2090::ItemOffset_t fixedSubzoneSize, bool                      sortItems); virtual ~OrbFESubzonePartitioner(); virtual ___465                  numCellsInZone() const; virtual ___2090::SubzoneOffset_t ___2783() const; virtual ___2090::ItemOffset_t    ___2782(___2090::SubzoneOffset_t ___469) const; virtual ___465                  ___4608(___2090 ___451) const; virtual ___2090                  szCoordinateAtZoneCell(___465 zoneCell) const; virtual ___2718                  numNodesInZone() const; virtual ___2090::SubzoneOffset_t ___2823() const; virtual ___2090::ItemOffset_t    ___2822(___2090::SubzoneOffset_t ___2734) const; virtual ___2718                  ___4657(___2090 nodeAddress) const; virtual ___2090                  ___3924(___2718 ___4656) const; virtual void                         setNodeSubzoneCoordinate(___2718 ___4656, ___2090 ___2759); }; }}