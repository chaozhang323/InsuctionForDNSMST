 #pragma once
#include "ThirdPartyHeadersBegin.h"
#   include <string>
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
#include "basicTypes.h"
#include "ClassMacros.h"
#include "FileIOStatistics.h"
#include "FileWriterInterface.h"
#include "mpiDatatype.h"
namespace tecplot { namespace teciompi { class MPIFileWriter : public ___3933::FileWriterInterface { public: MPIFileWriter( std::string const& ___1394, MPI_Comm           comm); virtual ~MPIFileWriter(); virtual ___372 ___2041() const; virtual ___372 open(); virtual ___372 close(bool ___3361); virtual ___3933::___1393 fileLoc(); virtual ___372 ___3460(); virtual ___372 ___3459(___3933::___1393 fileLoc); virtual std::string const& ___1394() const; virtual void ___3494(___372 ___2002); virtual ___372 ___2002() const; virtual void setDataFileType(DataFileType_e ___844); virtual DataFileType_e ___844() const; virtual ___3933::FileIOStatistics& statistics(); virtual size_t fwrite(void const* ___416, size_t size, size_t count); virtual int fprintf(char const* format, ...); private: ___3933::FileIOStatistics m_statistics; MPI_Comm                m_comm; MPI_File                m_fileHandle; std::string const       ___2461; bool                    m_isAscii; DataFileType_e          m_dataFileType; UNCOPYABLE_CLASS(MPIFileWriter); }; }}
