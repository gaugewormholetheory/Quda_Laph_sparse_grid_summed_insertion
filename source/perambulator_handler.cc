#include "perambulator_handler.h"
#include "multi_compare.h"
#include "stop_watch.h"
#include "field_ops.h"

#if defined(USE_GSL_CBLAS)
#include "gsl_cblas.h"
#elif defined(USE_OPENBLAS)
#include "cblas.h"
#endif

using namespace std;

typedef std::complex<double> dcmplx;
typedef std::complex<float>  fcmplx;

namespace LaphEnv {

PerambulatorHandler::PerambulatorType PerambulatorHandler::toPerambulatorType(const std::string& s)
   {
      if (s=="PERAMBULATOR_STANDARD"){
         printLaph(s + " output enabled."); return PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD;}                      
      if (s=="PERAMBULATOR_LEFT_SPARSE_GRID"){
         printLaph(s + " output enabled."); return PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID;}              
      if (s=="PERAMBULATOR_LEFT_INSERTION"){
         printLaph(s + " output enabled."); return PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION;}                  
      if (s=="PERAMBULATOR_RIGHT_INSERTION"){
         printLaph(s + " output enabled."); return PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION;}
      if (s=="PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION"){
         printLaph(s + " output enabled."); return PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION;} 
      if (s=="PERAMBULATOR_SUMMED_INSERTION") {
         printLaph(s + " output enabled."); return PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION;}                 
      if (s=="PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION"){
         printLaph(s + " output enabled."); return PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION;}
      return PerambulatorHandler::PerambulatorType::PERAMBULATOR_INVALID_ENUM;
   }
   
// *************************************************************************

void PerambulatorHandler::RecordKey::output(XMLHandler& xmlw) const 
{
 xmlw.set_root("RecordKey");
 xmlw.put_child("SinkSpin",make_string(getSinkSpin()));
 xmlw.put_child("SinkTime",make_string(getSinkTime()));
 xmlw.put_child("SourceLaphEigvecIndex",make_string(getSourceLaphEigvecIndex()));
}

std::string PerambulatorHandler::RecordKey::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}

// *************************************************************************


PerambulatorHandler::FileKey::FileKey(int in_srctime, int in_srcspin)
    : src_time(in_srctime), src_spin(in_srcspin) {}

PerambulatorHandler::FileKey::FileKey(XMLHandler& xmlr)
{
 try{
    XMLHandler xmlf(xmlr,"FileKey");
    xmlread(xmlf,"SourceTime",src_time,"PerambulatorHandler::FileKey");
    xmlread(xmlf,"SourceSpin",src_spin,"PerambulatorHandler::FileKey");}
 catch(...){
    errorLaph("Could not read PerambulatorHandler::FileKey");}
}

PerambulatorHandler::FileKey::FileKey(const FileKey& rhs)
    : src_time(rhs.src_time), src_spin(rhs.src_spin) {}

PerambulatorHandler::FileKey& PerambulatorHandler::FileKey::operator=(const FileKey& rhs)
{
 src_time=rhs.src_time;
 src_spin=rhs.src_spin;
 return *this;
}

void PerambulatorHandler::FileKey::output(XMLHandler& xmlw) const
{
 xmlw.set_root("FileKey");
 xmlw.put_child("SourceTime",make_string(src_time));
 xmlw.put_child("SourceSpin",make_string(src_spin));
}

std::string PerambulatorHandler::FileKey::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}

bool PerambulatorHandler::FileKey::operator<(const FileKey& rhs) const
{
 return multiLessThan(src_time,rhs.src_time,src_spin,rhs.src_spin);
}

bool PerambulatorHandler::FileKey::operator==(const FileKey& rhs) const
{
 return multiEqual(src_time,rhs.src_time,src_spin,rhs.src_spin);
}

bool PerambulatorHandler::FileKey::operator!=(const FileKey& rhs) const
{
 return multiNotEqual(src_time,rhs.src_time,src_spin,rhs.src_spin);
}


// *************************************************************************


// constructors: default, standard, sparse-grid-type, insertion-type, sparse-grid-insertion-type

PerambulatorHandler::PerambulatorHandler()
          : uPtr(0), gSmearPtr(0), qSmearPtr(0), qactionPtr(0),
            pmapPtr(0), invertPtr(0), Nspin(4), mode(ReadOnly),
            preconditioner(0), fPtr(0), fPtrSparseGrid(0), fPtrSummedInsertion(0), fPtrLeftInsertion(0), 
            fPtrRightInsertion(0), fPtrSparseGridSummedInsertion(0), fPtrSparseGridRightInsertion(0), 
            DHputPtr(0), DHputPtrSparseGrid(0), DHputPtrSummedInsertion(0), DHputPtrLeftInsertion(0), 
            DHputPtrRightInsertion(0), DHputPtrSparseGridSummedInsertion(0), DHputPtrSparseGridRightInsertion(0), 
            DHgetPtr(0), DHgetPtrSparseGrid(0), DHgetPtrSummedInsertion(0), DHgetPtrLeftInsertion(0), 
            DHgetPtrRightInsertion(0), DHgetPtrSparseGridSummedInsertion(0), DHgetPtrSparseGridRightInsertion(0) {}


PerambulatorHandler::PerambulatorHandler(const GaugeConfigurationInfo& gaugeinfo,
                                         const GluonSmearingInfo& gluonsmear,
                                         const QuarkSmearingInfo& quarksmear,
                                         const QuarkActionInfo& quark,
                                         const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                         const string& smeared_quark_filestub,
                                         bool upper_spin_components_only,
                                         Mode in_mode, const string& gauge_str)
          : invertPtr(0), preconditioner(0), DHputPtr(0), DHputPtrSparseGrid(0), DHputPtrSummedInsertion(0), 
            DHputPtrLeftInsertion(0), DHputPtrRightInsertion(0), DHputPtrSparseGridSummedInsertion(0), 
            DHputPtrSparseGridRightInsertion(0), DHgetPtr(0), DHgetPtrSparseGrid(0), DHgetPtrSummedInsertion(0), 
            DHgetPtrLeftInsertion(0), DHgetPtrRightInsertion(0), DHgetPtrSparseGridSummedInsertion(0), 
            DHgetPtrSparseGridRightInsertion(0){
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,pmap,smeared_quark_filestub,
		 upper_spin_components_only,gauge_str,in_mode);
}

PerambulatorHandler::PerambulatorHandler(const GaugeConfigurationInfo& gaugeinfo,
                                         const GluonSmearingInfo& gluonsmear,
                                         const QuarkSmearingInfo& quarksmear,
                                         const QuarkActionInfo& quark,
                                         const RandomSparseGrid& rsgrid,
                                         const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                         const string& smeared_quark_filestub,
                                         bool upper_spin_components_only,
                                         Mode in_mode, const string& gauge_str)
          : invertPtr(0), preconditioner(0), DHputPtr(0), DHputPtrSparseGrid(0), DHputPtrSummedInsertion(0), 
            DHputPtrLeftInsertion(0), DHputPtrRightInsertion(0), DHputPtrSparseGridSummedInsertion(0), 
            DHputPtrSparseGridRightInsertion(0), DHgetPtr(0), DHgetPtrSparseGrid(0), DHgetPtrSummedInsertion(0), 
            DHgetPtrLeftInsertion(0), DHgetPtrRightInsertion(0), DHgetPtrSparseGridSummedInsertion(0), 
            DHgetPtrSparseGridRightInsertion(0){
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,rsgrid,pmap,smeared_quark_filestub,
   upper_spin_components_only,gauge_str,in_mode);
}

PerambulatorHandler::PerambulatorHandler(const GaugeConfigurationInfo& gaugeinfo,
                                         const GluonSmearingInfo& gluonsmear,
                                         const QuarkSmearingInfo& quarksmear,
                                         const QuarkActionInfo& quark,
                                         const InsertionInfo& insinfo,
                                         const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                         const string& smeared_quark_filestub,
                                         bool upper_spin_components_only,
                                         Mode in_mode, const string& gauge_str)
          : invertPtr(0), preconditioner(0), DHputPtr(0), DHputPtrSparseGrid(0), DHputPtrSummedInsertion(0), 
            DHputPtrLeftInsertion(0), DHputPtrRightInsertion(0), DHputPtrSparseGridSummedInsertion(0), 
            DHputPtrSparseGridRightInsertion(0), DHgetPtr(0), DHgetPtrSparseGrid(0), DHgetPtrSummedInsertion(0), 
            DHgetPtrLeftInsertion(0), DHgetPtrRightInsertion(0), DHgetPtrSparseGridSummedInsertion(0), 
            DHgetPtrSparseGridRightInsertion(0){
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,insinfo,pmap,smeared_quark_filestub,
   upper_spin_components_only,gauge_str,in_mode);
}

PerambulatorHandler::PerambulatorHandler(const GaugeConfigurationInfo& gaugeinfo,
                                         const GluonSmearingInfo& gluonsmear,
                                         const QuarkSmearingInfo& quarksmear,
                                         const QuarkActionInfo& quark,
                                         const RandomSparseGrid& rsgrid,
                                         const InsertionInfo& insinfo,
                                         const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                         const string& smeared_quark_filestub,
                                         bool upper_spin_components_only,
                                         Mode in_mode, const string& gauge_str)
          : invertPtr(0), preconditioner(0), DHputPtr(0), DHputPtrSparseGrid(0), DHputPtrSummedInsertion(0), 
            DHputPtrLeftInsertion(0), DHputPtrRightInsertion(0), DHputPtrSparseGridSummedInsertion(0), 
            DHputPtrSparseGridRightInsertion(0), DHgetPtr(0), DHgetPtrSparseGrid(0), DHgetPtrSummedInsertion(0), 
            DHgetPtrLeftInsertion(0), DHgetPtrRightInsertion(0), DHgetPtrSparseGridSummedInsertion(0), 
            DHgetPtrSparseGridRightInsertion(0){
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,rsgrid,insinfo,pmap,
          smeared_quark_filestub,upper_spin_components_only,gauge_str,in_mode);
}

void PerambulatorHandler::setInfo(const GaugeConfigurationInfo& gaugeinfo,
                                  const GluonSmearingInfo& gluonsmear,
                                  const QuarkSmearingInfo& quarksmear,
                                  const QuarkActionInfo& quark,
				  const RandomSparseGrid& rsgrid,
                                  const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                  const string& smeared_quark_filestub,
                                  bool upper_spin_components_only,
                                  Mode in_mode, const string& gauge_str)
{
 clear(); 
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,rsgrid,pmap,
          smeared_quark_filestub,upper_spin_components_only,gauge_str,in_mode);
}

void PerambulatorHandler::setInfo(const GaugeConfigurationInfo& gaugeinfo,
                                  const GluonSmearingInfo& gluonsmear,
                                  const QuarkSmearingInfo& quarksmear,
                                  const QuarkActionInfo& quark,
                                  const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                  const string& smeared_quark_filestub,
                                  bool upper_spin_components_only,
                                  Mode in_mode, const string& gauge_str)
{
 clear(); 
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,pmap,smeared_quark_filestub,
		 upper_spin_components_only,gauge_str,in_mode);
}

void PerambulatorHandler::setInfo(const GaugeConfigurationInfo& gaugeinfo,
                                  const GluonSmearingInfo& gluonsmear,
                                  const QuarkSmearingInfo& quarksmear,
                                  const QuarkActionInfo& quark,
			  	  const InsertionInfo& insinfo,
                                  const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                  const string& smeared_quark_filestub,
                                  bool upper_spin_components_only,
                                  Mode in_mode, const string& gauge_str)
{
 clear(); 
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,insinfo,pmap,
          smeared_quark_filestub,upper_spin_components_only,gauge_str,in_mode);
}

void PerambulatorHandler::setInfo(const GaugeConfigurationInfo& gaugeinfo,
                                  const GluonSmearingInfo& gluonsmear,
                                  const QuarkSmearingInfo& quarksmear,
                                  const QuarkActionInfo& quark,
                                  const RandomSparseGrid& rsgrid,
											 const InsertionInfo& insinfo,
                                  const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                  const string& smeared_quark_filestub,
                                  bool upper_spin_components_only,
                                  Mode in_mode, const string& gauge_str)
{
 clear(); 
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,rsgrid,insinfo,pmap,
          smeared_quark_filestub,upper_spin_components_only,gauge_str,in_mode);
}

void PerambulatorHandler::set_info(const GaugeConfigurationInfo& gaugeinfo,
                                   const GluonSmearingInfo& gluonsmear,
                                   const QuarkSmearingInfo& quarksmear,
                                   const QuarkActionInfo& quark,
											  const RandomSparseGrid& rsgrid,
                                   const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                   const string& smeared_quark_filestub,
                                   bool upper_spin_components_only,
                                   const string& gauge_str, Mode in_mode)
{
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,pmap,smeared_quark_filestub,
		 upper_spin_components_only,gauge_str,in_mode);
 try{

	 sgHandler = new SparseGridHandler(gaugeinfo,gluonsmear,quarksmear,quark,rsgrid,Nspin);
    if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID)==1){
       fPtrSparseGrid = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID) );
       if ((mode==Compute)||(mode==Merge)){
          DHputPtrSparseGrid=new DataPutHandlerMF<SparseGridHandler,FileKey,RecordKey,DataType>(
                       *sgHandler,*fPtrSparseGrid,"Laph--SparseGridQuarkPeramb","PerambulatorHandlerSparseGridDataFile");
	   	}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtrSparseGrid=new DataGetHandlerMF<SparseGridHandler,FileKey,RecordKey,DataType>(
                       *sgHandler,*fPtrSparseGrid,"Laph--SparseGridQuarkPeramb","PerambulatorHandlerSparseGridDataFile",globalmode);
	    	}
      }
    }
 catch(const std::exception& xp){
    errorLaph(make_strf("allocation problem in PerambulatorHandler: %s",xp.what()));}
}

void PerambulatorHandler::set_info(const GaugeConfigurationInfo& gaugeinfo,
                                   const GluonSmearingInfo& gluonsmear,
                                   const QuarkSmearingInfo& quarksmear,
                                   const QuarkActionInfo& quark,
				   const InsertionInfo& insinfo,
                                   const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                   const string& smeared_quark_filestub,
                                   bool upper_spin_components_only,
                                   const string& gauge_str, Mode in_mode)
{
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,pmap,smeared_quark_filestub,
		 upper_spin_components_only,gauge_str,in_mode);

 try{    
	 insHandler = new InsertionHandler(gaugeinfo,gluonsmear,quarksmear,quark,insinfo,Nspin);
    sginsHandler = nullptr;
	 if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)==1){
       fPtrSummedInsertion = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION) );
	    if ((mode==Compute)||(mode==Merge)){
			 DHputPtrSummedInsertion=new DataPutHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrSummedInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile");}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtrSummedInsertion=new DataGetHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrSummedInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile",globalmode);}
       }  
	 if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION)==1){
       fPtrLeftInsertion = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION) );
	    if ((mode==Compute)||(mode==Merge)){
			 DHputPtrLeftInsertion=new DataPutHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrLeftInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile");}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtrLeftInsertion=new DataGetHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrLeftInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile",globalmode);}
       }
	 if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)==1){
       fPtrRightInsertion = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION) );
	    if ((mode==Compute)||(mode==Merge)){
			 DHputPtrRightInsertion=new DataPutHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrRightInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile");}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtrRightInsertion=new DataGetHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrRightInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile",globalmode);}
       }
 }
 catch(const std::exception& xp){
    errorLaph(make_strf("allocation problem in PerambulatorHandler: %s",xp.what()));}
}

void PerambulatorHandler::set_info(const GaugeConfigurationInfo& gaugeinfo,
                                   const GluonSmearingInfo& gluonsmear,
                                   const QuarkSmearingInfo& quarksmear,
                                   const QuarkActionInfo& quark,
											  const RandomSparseGrid& rsgrid,
                                   const InsertionInfo& insinfo,
                                   const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                   const string& smeared_quark_filestub,
                                   bool upper_spin_components_only,
                                   const string& gauge_str, Mode in_mode)
{
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,pmap,smeared_quark_filestub,
		 upper_spin_components_only,gauge_str,in_mode);
 try{

	 sgHandler = new SparseGridHandler(gaugeinfo,gluonsmear,quarksmear,quark,rsgrid,Nspin);
    insHandler = new InsertionHandler(gaugeinfo,gluonsmear,quarksmear,quark,insinfo,Nspin);
    sginsHandler = new SparseGridInsertionHandler(gaugeinfo,gluonsmear,quarksmear,quark,rsgrid,insinfo,Nspin);
    if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID)==1){
       fPtrSparseGrid = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID) );
	    if ((mode==Compute)||(mode==Merge)){
			 DHputPtrSparseGrid=new DataPutHandlerMF<SparseGridHandler,FileKey,RecordKey,DataType>(
                       *sgHandler,*fPtrSparseGrid,"Laph--SparseGridQuarkPeramb","PerambulatorHandlerSparseGridDataFile");}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtrSparseGrid=new DataGetHandlerMF<SparseGridHandler,FileKey,RecordKey,DataType>(
                       *sgHandler,*fPtrSparseGrid,"Laph--SparseGridQuarkPeramb","PerambulatorHandlerSparseGridDataFile",globalmode);}
       }  
	 if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)==1){
       fPtrSummedInsertion = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION) );
	    if ((mode==Compute)||(mode==Merge)){
			 DHputPtrSummedInsertion=new DataPutHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrSummedInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile");}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtrSummedInsertion=new DataGetHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrSummedInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile",globalmode);}
       }  
	 if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION)==1){
       fPtrLeftInsertion = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION) );
	    if ((mode==Compute)||(mode==Merge)){
			 DHputPtrLeftInsertion=new DataPutHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrLeftInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile");}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtrLeftInsertion=new DataGetHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrLeftInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile",globalmode);}
       }
	 if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)==1){
       fPtrRightInsertion = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION) );
	    if ((mode==Compute)||(mode==Merge)){
			 DHputPtrRightInsertion=new DataPutHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrRightInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile");}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtrRightInsertion=new DataGetHandlerMF<InsertionHandler,FileKey,RecordKey,DataType>(
                       *insHandler,*fPtrRightInsertion,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile",globalmode);}
       }
	 if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION)==1){
       fPtrSparseGridSummedInsertion = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION) );
	    if ((mode==Compute)||(mode==Merge)){
			 DHputPtrSparseGridSummedInsertion=new DataPutHandlerMF<SparseGridInsertionHandler,FileKey,RecordKey,DataType>(
                       *sginsHandler,*fPtrSparseGridSummedInsertion,"Laph--SparseGridQuarkPeramb","PerambulatorHandlerSparseGridDataFile");}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtrSparseGridSummedInsertion=new DataGetHandlerMF<SparseGridInsertionHandler,FileKey,RecordKey,DataType>(
                       *sginsHandler,*fPtrSparseGridSummedInsertion,"Laph--SparseGridQuarkPeramb","PerambulatorHandlerSparseGridDataFile",globalmode);}
       }
	 if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION)==1){
       fPtrSparseGridRightInsertion = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION) );
	    if ((mode==Compute)||(mode==Merge)){
			 DHputPtrSparseGridRightInsertion=new DataPutHandlerMF<SparseGridInsertionHandler,FileKey,RecordKey,DataType>(
                       *sginsHandler,*fPtrSparseGridRightInsertion,"Laph--SparseGridQuarkPeramb","PerambulatorHandlerSparseGridDataFile");}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtrSparseGridRightInsertion=new DataGetHandlerMF<SparseGridInsertionHandler,FileKey,RecordKey,DataType>(
                       *sginsHandler,*fPtrSparseGridRightInsertion,"Laph--SparseGridQuarkPeramb","PerambulatorHandlerSparseGridDataFile",globalmode);}
       }
 }
 catch(const std::exception& xp){
    errorLaph(make_strf("allocation problem in PerambulatorHandler: %s",xp.what()));}
}

void PerambulatorHandler::set_info(const GaugeConfigurationInfo& gaugeinfo,
                                   const GluonSmearingInfo& gluonsmear,
                                   const QuarkSmearingInfo& quarksmear,
                                   const QuarkActionInfo& quark,
                                   const std::map<PerambulatorHandler::PerambulatorType, FileListInfo>& pmap,
                                   const string& smeared_quark_filestub,
                                   bool upper_spin_components_only,
                                   const string& gauge_str, Mode in_mode)
{
 try{
    uPtr = new GaugeConfigurationInfo(gaugeinfo);
    gSmearPtr = new GluonSmearingInfo(gluonsmear);
    qSmearPtr = new QuarkSmearingInfo(quarksmear);
    qactionPtr = new QuarkActionInfo(quark);
	 Nspin = (upper_spin_components_only) ? 2 : 4;
    mode = in_mode;
    pmapPtr = new std::map<PerambulatorType, FileListInfo>(pmap);
    if ( pmap.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD)==1){
       fPtr = new FileListInfo( pmap.at(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD) );
       if ((mode==Compute)||(mode==Merge)){
          DHputPtr=new DataPutHandlerMF<PerambulatorHandler,FileKey,RecordKey,DataType>(
                          *this,*fPtr,"Laph--QuarkPeramb","PerambulatorHandlerDataFile");}
       if ((mode==ReadOnly)||(mode==Check)){
          bool globalmode=(mode==ReadOnly); 
          DHgetPtr=new DataGetHandlerMF<PerambulatorHandler,FileKey,RecordKey,DataType>(
                    *this,*fPtr,"Laph--QuarkPeramb","PerambulatorHandlerDataFile",globalmode);}
       }
    }  
 catch(const std::exception& xp){
    errorLaph(make_strf("allocation problem in PerambulatorHandler: %s",xp.what()));}

 if (mode==Compute){
    connectGaugeConfigurationHandler();
    connectQuarkSmearingHandler(smeared_quark_filestub);}
}

PerambulatorHandler::~PerambulatorHandler()
{
 clear();
}


void PerambulatorHandler::clear()
{
	if (isInfoSet()) {
		try{
			delete uPtr;
			delete gSmearPtr;
			delete qSmearPtr;
			delete qactionPtr;
			if (QudaInfo::clover_on_device){
				freeCloverQuda();
				QudaInfo::clover_on_device=false;}  // delete clover in case next task uses different kappa
			if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD) == 1){
            delete fPtr; fPtr=0;
            delete DHgetPtr; DHgetPtr=0;
		      delete DHputPtr; DHputPtr=0;}
			if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID) == 1){
            delete fPtrSparseGrid; fPtrSparseGrid=0;
            delete DHgetPtrSparseGrid; DHgetPtrSparseGrid=0;
		      delete DHputPtrSparseGrid; DHputPtrSparseGrid=0;}	
         if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION) == 1){
            delete fPtrSummedInsertion; fPtrSummedInsertion=0;
            delete DHgetPtrSummedInsertion; DHgetPtrSummedInsertion=0;
		      delete DHputPtrSummedInsertion; DHputPtrSummedInsertion=0;}
         if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION) == 1){
            delete fPtrLeftInsertion; fPtrLeftInsertion=0;
            delete DHgetPtrLeftInsertion; DHgetPtrLeftInsertion=0;
		      delete DHputPtrLeftInsertion; DHputPtrLeftInsertion=0;}
         if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION) == 1){
            delete fPtrRightInsertion; fPtrRightInsertion=0;
            delete DHgetPtrRightInsertion; DHgetPtrRightInsertion=0;
		      delete DHputPtrRightInsertion; DHputPtrRightInsertion=0;}
         if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION) == 1){
            delete fPtrSparseGridSummedInsertion; fPtrSparseGridSummedInsertion=0;
            delete DHgetPtrSparseGridSummedInsertion; DHgetPtrSparseGridSummedInsertion=0;
		      delete DHputPtrSparseGridSummedInsertion; DHputPtrSparseGridSummedInsertion=0;}
         if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION) == 1){
            delete fPtrSparseGridRightInsertion; fPtrSparseGridRightInsertion=0;
            delete DHgetPtrSparseGridRightInsertion; DHgetPtrSparseGridRightInsertion=0;
		      delete DHputPtrSparseGridRightInsertion; DHputPtrSparseGridRightInsertion=0;}
         delete pmapPtr; delete pmapPtr;
			delete invertPtr;
			if (preconditioner){
				destroyMultigridQuda(preconditioner);}
			clearComputationSet();}
		catch(const std::exception& xp){
			errorLaph("abort");}
		uPtr=0;
		gSmearPtr=0;
		qSmearPtr=0;
		qactionPtr=0;
		invertPtr=0;
		mode=ReadOnly;
		preconditioner=0;
		disconnectGaugeConfigurationHandler();
		disconnectQuarkSmearingHandler();

	}
}


  // ********************************
  // *
  // *    sub-handler connections  (private)
  // *
  // ********************************

void PerambulatorHandler::connectGaugeConfigurationHandler()
{
 if (gaugeHandler.get()==0){
    try{
       gaugeHandler.reset(new GaugeConfigurationHandler(*uPtr));}
    catch(const std::exception& xp){
       errorLaph("allocation problem in PerambulatorHandler::connectGaugeConfigurationHandler");}}
}

void PerambulatorHandler::disconnectGaugeConfigurationHandler()
{
 try{ gaugeHandler.reset();}
 catch(const std::exception& xp){
    errorLaph("delete problem in PerambulatorHandler::disconnectGluonSmearingHandler");}
}

void PerambulatorHandler::connectQuarkSmearingHandler(const string& smeared_quark_filestub)
{
 if (qSmearHandler.get()==0){
    try{
       qSmearHandler.reset(new QuarkSmearingHandler(*gSmearPtr,*uPtr,*qSmearPtr,
                                                    smeared_quark_filestub));}
    catch(const std::exception& xp){
       errorLaph("allocation problem in PerambulatorHandler::connectQuarkSmearingHandler");}}
}

void PerambulatorHandler::disconnectQuarkSmearingHandler()
{
 try{ qSmearHandler.reset();}
 catch(const std::exception& xp){
    errorLaph("delete problem in PerambulatorHandler::disconnectQuarkSmearingHandler");}
}



void PerambulatorHandler::setInverter(const InverterInfo& invinfo)
{
 if (mode!=Compute){
    errorLaph("Cannot setInverter in PerambulatorHandler unless in compute mode");}
 if (!isInfoSet()){
    errorLaph("Cannot setInverter in PerambulatorHandler unless info is set");}
 try{
    delete invertPtr;
    invertPtr = new InverterInfo(invinfo);
    }
 catch(const std::exception& xp){
    errorLaph("allocation error in PerambulatorHandler::setInverter");}

     // set up the inverter params for quda
 invertPtr->setQudaInvertParam(quda_inv_param,*qactionPtr);

     // check that Dirac-Pauli basis is used
 if (quda_inv_param.gamma_basis!=QUDA_DIRAC_PAULI_GAMMA_BASIS){
    errorQuda("The Dirac-Pauli basis must be used");}
}

const InverterInfo& PerambulatorHandler::getInverterInfo() const 
{
 if (invertPtr!=0){
    errorLaph("error in PerambulatorHandler: must setInverter before calling getInverterInfo");}
 return *invertPtr;
}


void PerambulatorHandler::clearComputationSet()
{
 perambComps.computations.clear();
}


void PerambulatorHandler::setComputationSet(const XMLHandler& xmlin)
{
 if (mode!=Compute){
    errorLaph("Cannot setComputationSet in PerambulatorHandler unless in compute mode");}
 if (!isInfoSet()){
    errorLaph("Cannot setComputationSet in PerambulatorHandler unless info is set");}
 XMLHandler xmlrdr(xmlin);
 if (!perambComps.computations.empty()) perambComps.computations.clear();
 uint Textent = uPtr->getTimeExtent();
 uint nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 if (xml_tag_count(xmlrdr,"ComputationSet")!=1){
    errorLaph("Cannot setComputationSet in PerambulatorHandler since no <ComputationSet> tag");}
 
 XMLHandler xmlrd(xmlrdr,"ComputationSet");
 string multisrc_reply="yes";
 xmlreadif(xmlrd,"UseMultiSrcInverter",multisrc_reply,"LAPH_PERAMBULATORS");
 if ((multisrc_reply=="YES")||(multisrc_reply=="Y")||(multisrc_reply=="Yes")||(multisrc_reply=="y")){
    multisrc_reply="yes";}
 if (multisrc_reply!="yes") multisrc_reply="no";
 uint nSinkLaphBatch,nSinkQudaBatch,nEigQudaBatch,nSinkTimes;
 xmlread(xmlrd,"NumSinksBeforeProject",nSinkLaphBatch,"LAPH_PERAMBULATORS");
 xmlread(xmlrd,"NumSinksInProjectBatch",nSinkQudaBatch,"LAPH_PERAMBULATORS"); 
 xmlread(xmlrd,"NumEigsInProjectBatch",nEigQudaBatch,"LAPH_PERAMBULATORS");

 if ((nSinkLaphBatch > nEigs*Nspin) || (nSinkLaphBatch == 0)){
    errorLaph(make_strf("Invalid value %d for <NumSinksBeforeProject>: must be between 1 and nEigs*Nspin=%d",
              nSinkLaphBatch,nEigs*Nspin));}
 if ((nEigQudaBatch > nEigs) || (nEigQudaBatch == 0)){
    errorLaph(make_strf("Invalid value %d for <NumEigsInProjectBatch>: must be between 1 and nEigs=%d",
              nEigQudaBatch,nEigs));}
 if ((nSinkQudaBatch > nSinkLaphBatch) || (nSinkQudaBatch == 0)){
    errorLaph(make_strf("Invalid value %d for <NumSinksInProjectBatch>: must be between 1 and %d",
              nSinkQudaBatch,nSinkLaphBatch));}
 perambComps.nSinkLaphBatch=nSinkLaphBatch;
 perambComps.nSinkQudaBatch=nSinkQudaBatch;
 perambComps.nEigQudaBatch=nEigQudaBatch;
 perambComps.useMultiSrcInverter=(multisrc_reply=="yes")?true:false;

 list<XMLHandler> xmlcs(xmlrd.find("Computation"));
 for (list<XMLHandler>::iterator it=xmlcs.begin();it!=xmlcs.end();++it){
    int source_time, tmin, tmax;
    xmlread(*it,"SourceTime",source_time,"LAPH_PERAMBULATORS");
    tmin=0; tmax=Textent-1;  
		xmlreadif(*it,"SinkMinimumTime",tmin,"LAPH_PERAMBULATORS");
		xmlreadif(*it,"SinkMaximumTime",tmax,"LAPH_PERAMBULATORS");
    if ((source_time<0)||(source_time>=int(Textent))){
       errorLaph(make_strf("Invalid source time %d",source_time));}
		int nSinkTimes=tmax-tmin+1; 
    if ((tmin<=(-int(Textent)))||(tmax>=int(Textent))||(nSinkTimes<=0)||(nSinkTimes>int(Textent))){
       errorLaph(make_strf("Invalid sink tmin and tmax (%d,%d)",tmin,tmax));}
    set<int> srcev_indices;
    if (xml_tag_count(*it,"SourceLaphEigvecIndices")==1){
       vector<int> srcev_inds;
       xmlread(*it,"SourceLaphEigvecIndices",srcev_inds,"LAPH_PERAMBULATORS");
       for (int k=0;k<int(srcev_inds.size());++k){
          if ((srcev_inds[k]<0)||(srcev_inds[k]>=int(nEigs))){
             errorLaph(make_strf("Invalid src laph eigvec index %d",srcev_inds[k]));}
       srcev_indices.insert(srcev_inds[k]);}}
    else if (  (xml_tag_count(*it,"SourceLaphEigvecIndexMin")==1)
             &&(xml_tag_count(*it,"SourceLaphEigvecIndexMax")==1)){
       int sevmin=-1, sevmax=-1;
       xmlread(*it,"SourceLaphEigvecIndexMin",sevmin,"LAPH_PERAMBULATORS");
       xmlread(*it,"SourceLaphEigvecIndexMax",sevmax,"LAPH_PERAMBULATORS");
       if (sevmin<0) sevmin=0;
       if (sevmax>=int(nEigs)) sevmax=nEigs-1;
       if (sevmax<sevmin) sevmax=sevmin;
       for (int ev=sevmin;ev<=sevmax;++ev){
          srcev_indices.insert(ev);}}
    else{
       for (int ev=0;ev<int(nEigs);++ev){
          srcev_indices.insert(ev);}}
    if (srcev_indices.empty()){
       errorLaph("Empty src laph eigvec indices set");}
    perambComps.computations.push_back(
          PerambComputation(source_time,tmin,tmax,srcev_indices));}

 printLaph("\nLAPH_PERAMBULATORS sink computations:\n");
 printLaph(make_strf("  NumSinksBeforeProject = %d",perambComps.nSinkLaphBatch));
 printLaph(make_strf(" NumSinksInProjectBatch = %d",perambComps.nSinkQudaBatch));
 printLaph(make_strf("  NumEigsInProjectBatch = %d\n",perambComps.nEigQudaBatch));
 printLaph(make_strf("    UseMultiSrcInverter = %s\n",multisrc_reply));
 printLaph(make_strf(" Number of computations = %d",perambComps.computations.size()));
 int count=0;
 for (list<PerambComputation>::const_iterator it=perambComps.computations.begin();
      it!=perambComps.computations.end();count++,it++){
    XMLHandler xmlout("Computation");
    xmlout.put_child("SourceTime",make_string(it->src_time));
    xmlout.put_child("SinkMinimumTime",make_string(it->tmin));
    xmlout.put_child("SinkMaximumTime",make_string(it->tmax));
    string srcevindstr;
    const set<int>& indices(it->src_lapheigvec_indices);
    set<int>::const_iterator vt=indices.begin();
    int rangestart=*vt; int rangestop=*vt+1; ++vt;
    bool keep_going=true;
    while (keep_going){
       if ((vt!=indices.end())&&(rangestop==(*vt))){   // in a started range, increase rangestop
          rangestop++; ++vt;}
       else{                    // end of a range
          if (rangestop==(rangestart+1)){             // only one in range
             srcevindstr+=make_string(rangestart);}
          else{                                       // more than one in range
             srcevindstr+=make_string(rangestart)+"-"+make_string(rangestop-1);}
          if (vt!=indices.end()){
             rangestart=*vt; rangestop=rangestart+1;
             srcevindstr+=" "; ++vt;}
          else{
             keep_going=false;}}}
    xmlout.put_child("SourceLaphEigvecIndices",srcevindstr);
    printLaph(make_strf("\nComputation %d:",count));
    printLaph(make_strf("%s",xmlout.output()));}
 printLaph("\n");
}


void PerambulatorHandler::getFileMap(XMLHandler& xmlout) const
{
 if (isInfoSet()) DHputPtr->getFileMap(xmlout);
}

map<int,PerambulatorHandler::FileKey> PerambulatorHandler::getSuffixMap() const
{
 check_info_set("getSuffixMap");
 if (mode!=ReadOnly){
    return DHputPtr->getSuffixMap();}
 else{
    return DHgetPtr->getSuffixMap();}
}


void PerambulatorHandler::outputSuffixMap()
{
 check_info_set("getSuffixMap");
 map<int,PerambulatorHandler::FileKey> suffixmap=getSuffixMap();
 printLaph("\nSuffix map:");
 for (map<int,PerambulatorHandler::FileKey>::const_iterator it=suffixmap.begin();
      it!=suffixmap.end();++it){
    printLaph(make_strf("suffix  %d:  source time = %d  source spin = %d",
               it->first,it->second.src_time,it->second.src_spin));}
 printLaph("\n");
}


void PerambulatorHandler::outputSuffixMap(XMLHandler& xmlout)
{
 check_info_set("getSuffixMap");
 map<int,PerambulatorHandler::FileKey> suffixmap=getSuffixMap();
 xmlout.set_root("SuffixMap");
 for (map<int,PerambulatorHandler::FileKey>::const_iterator it=suffixmap.begin();
      it!=suffixmap.end();++it){
    XMLHandler xmltmp;
    xmltmp.set_root("Suffix");
    xmltmp.put_child("Index",make_string(it->first));
    xmltmp.put_child("SourceTime",make_string(it->second.src_time));
    xmltmp.put_child("SourceSpin",make_string(it->second.src_spin));
    xmlout.put_child(xmltmp);}
}


bool PerambulatorHandler::isInfoSet() const
{
 return ((uPtr!=0)&&(gSmearPtr!=0)&&(qSmearPtr!=0)//&&(fPtr!=0)
        &&(qactionPtr!=0));
}


void PerambulatorHandler::check_info_set(const string& name) const
{
 if (!isInfoSet()){
    errorLaph(make_strf("error in PerambulatorHandler: must setInfo before calling %s",name));}
}


const GaugeConfigurationInfo& PerambulatorHandler::getGaugeConfigurationInfo() const 
{
 check_info_set("getGaugeConfigurationInfo");
 return *uPtr;
}

const GluonSmearingInfo& PerambulatorHandler::getGluonSmearingInfo() const
{
 check_info_set("getGluonSmearingInfo");
 return *gSmearPtr;
}

const QuarkSmearingInfo& PerambulatorHandler::getQuarkSmearingInfo() const
{
 check_info_set("getQuarkSmearingInfo");
 return *qSmearPtr;
}

const QuarkActionInfo& PerambulatorHandler::getQuarkActionInfo() const 
{
 check_info_set("getQuarkActionInfo");
 return *qactionPtr;
}

const FileListInfo& PerambulatorHandler::getFileListInfo() const 
{ if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD)!=1)
	errorLaph("There is no standard file list.");
 check_info_set("getFileListInfo");
 return *fPtr;
}

const FileListInfo& PerambulatorHandler::getFileListInfoSparseGrid() const 
{
 if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID)!=1)
    errorLaph("There is no sparse grid file list.");
 check_info_set("getFileListInfoSparseGrid"); 
 return *fPtrSparseGrid;
}

const FileListInfo& PerambulatorHandler::getFileListInfoSummedInsertion() const 
{
 if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)!=1)
    errorLaph("There is no summed insertion file list.");
 check_info_set("getFileListInfoSummedInsertion");
 return *fPtrSummedInsertion;
}

const FileListInfo& PerambulatorHandler::getFileListInfoLeftInsertion() const
{
 if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION)!=1)
     errorLaph("There is no left insertion file list.");
 check_info_set("getFileListInfoLeftInsertion");
 return *fPtrLeftInsertion;
}

const FileListInfo& PerambulatorHandler::getFileListInfoRightInsertion() const
{
 if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)!=1)
     errorLaph("There is no right insertion file list.");
 check_info_set("getFileListInfoRightInsertion");
 return *fPtrRightInsertion;
}

const FileListInfo& PerambulatorHandler::getFileListInfoSparseGridSummedInsertion() const 
{
if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION)!=1)    
     errorLaph("There is no sparse grid summed insertion file list.");
 check_info_set("getFileListInfoSparseGridSummedInsertion");
 return *fPtrSparseGridSummedInsertion;
}

const FileListInfo& PerambulatorHandler::getFileListInfoSparseGridRightInsertion() const
{       
 if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION)!=1)
      errorLaph("There is no sparse grid right insertion file list.");
 check_info_set("getFileListInfoSparseGridRightInsertion");
 return *fPtrSparseGridRightInsertion;
}

bool PerambulatorHandler::standardOutput() const {
        return ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD)!=1);         
}

bool PerambulatorHandler::sparseGridOutput() const { 
	return ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID)!=1); 
}

bool PerambulatorHandler::summedInsertionOutput() const {
   return ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)!=1);
}

bool PerambulatorHandler::leftInsertionOutput() const {
   return ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION)!=1);
}

bool PerambulatorHandler::rightInsertionOutput() const {
   return ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)!=1);
}

bool PerambulatorHandler::sparseGridsummedInsertionOutput() const {
   return ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION)!=1);
}

bool PerambulatorHandler::sparseGridrightInsertionOutput() const {
    return ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION)!=1);
}

uint PerambulatorHandler::getNumberOfLaplacianEigenvectors() const
{
 return qSmearPtr->getNumberOfLaplacianEigenvectors();
}

int PerambulatorHandler::getTimeExtent() const 
{
 check_info_set("getTimeExtent");
 return uPtr->getTimeExtent();
}

int PerambulatorHandler::getNSpin() const 
{
 check_info_set("getNSpin");
 return Nspin;
}

void PerambulatorHandler::getHeader(XMLHandler& xmlout) const
{
 check_info_set("getHeader"); 
 xmlout.set_root("PerambulatorHandlerDataFile");
 XMLHandler xmltmp;
 uPtr->output(xmltmp); xmlout.put_child(xmltmp);
 gSmearPtr->output(xmltmp); xmlout.put_child(xmltmp);
 qSmearPtr->output(xmltmp); xmlout.put_child(xmltmp);
 qactionPtr->output(xmltmp); xmlout.put_child(xmltmp);
 xmlout.put_child("NumSpinComponents",make_string(Nspin));
}


bool PerambulatorHandler::checkHeader(XMLHandler& xmlin, int suffix)
{
 check_info_set("getHeader"); 
 XMLHandler xml_in(xmlin);
 if (xml_tag_count(xml_in,"PerambulatorHandlerDataFile")!=1) return false;
 XMLHandler xmlr(xml_in,"PerambulatorHandlerDataFile");
 GaugeConfigurationInfo gauge_check(xmlr);
 GluonSmearingInfo gsmear_check(xmlr);
 QuarkSmearingInfo qsmear_check(xmlr);
 uint numspin;
 xmlread(xmlr,"NumSpinComponents", numspin, "PerambulatorHandler");
 QuarkActionInfo qaction_check(xmlr);
 try {
    uPtr->checkEqual(gauge_check);
    gSmearPtr->checkEqual(gsmear_check);
    qSmearPtr->checkEqual(qsmear_check); 
    if (numspin!=Nspin){
       throw(std::invalid_argument("Perambulator checkEqual failed...NumSpinComponents mismatch"));}
    qactionPtr->checkEqual(qaction_check); }
 catch(const exception& xp){ return false;}
 return true;
}

void PerambulatorHandler::writeHeader(XMLHandler& xmlout, 
                                     const PerambulatorHandler::FileKey& fkey,
                                     int suffix)
{
 check_info_set("writeHeader"); 
 xmlout.set_root("PerambulatorHandlerDataFile");
 XMLHandler xmltmp;
 uPtr->output(xmltmp); xmlout.put_child(xmltmp);
 gSmearPtr->output(xmltmp); xmlout.put_child(xmltmp);
 qSmearPtr->output(xmltmp); xmlout.put_child(xmltmp);
 qactionPtr->output(xmltmp); xmlout.put_child(xmltmp);
 fkey.output(xmltmp); xmlout.put_child(xmltmp);
 xmlout.put_child("NumSpinComponents",make_string(Nspin));
}

 // ************************************************************************************


bool PerambulatorHandler::setUpPreconditioning(QudaInvertParam& invParam)
{
 if (mode!=Compute){
    errorLaph("Cannot setUpPreconditioning in PerambulatorHandler unless in compute mode");}
 if (!isInfoSet()){
    errorLaph("Cannot setUpPreconditioning in PerambulatorHandler unless info is set");}
 if (invertPtr->getName()=="GCR_MULTIGRID"){
    if (invertPtr->QudaMGInfoPtr){
       if (preconditioner){
          destroyMultigridQuda(preconditioner);}
       printLaph("Initializing the fabulous Multigrid preconditioner!");
       preconditioner=newMultigridQuda(&invertPtr->QudaMGInfoPtr->mg_param);
       invParam.preconditioner=preconditioner;
       return true;}
    else{
       throw(std::invalid_argument("setUpPreconditioning failed since invParam not set"));}}
 return false;
}


 // ************************************************************************************

    //  do "perambulator" inversions for all spin indices, one time source and a set of source 
    //  laph eigenvector numbers.  All sink spin indices, sink times, and laph_eigvec indices 
    //  at the sink are computed. 
 

void PerambulatorHandler::computePerambulators(bool extra_soln_check, bool print_coeffs, bool report_gflops)
{
 if ((!isInfoSet())||(invertPtr==0)||(mode!=Compute)){
    errorLaph(make_str("cannot computePerambulators in PerambulatorHandler until info ",
              "and inverter set and in compute mode"));}
 StopWatch totaltime; totaltime.start();

     // check temporal boundary conditions in InverterInfo and GaugeConfigurationInfo
     // are the same
 bool tbc1=qactionPtr->isFermionTimeBCAntiPeriodic();
 bool tbc2=uPtr->isFermionTimeBCAntiPeriodic();
 if (tbc1!=tbc2){
    errorLaph("Inconsistent fermion time boundary conditions in QuarkActionInfo and GaugeConfigurationInfo",true);}

     // load the gauge configuration onto host and device
 StopWatch rolex; rolex.start();
 gaugeHandler->setData();
 XMLHandler gauge_xmlinfo;
 gaugeHandler->getXMLInfo(gauge_xmlinfo);
 printLaph("XML info for the gauge configuration:");
 printLaph(make_strf("%s\n",gauge_xmlinfo.output()));
 gaugeHandler->copyDataToDevice();
 rolex.stop();
 double grtime=rolex.getTimeInSeconds();
 printLaph("...gauge configuration loaded on host and device");
 printLaph(make_str(" time to load gauge configuration was ",grtime," seconds"));

     // load the clover term if needed
 double clovertime=0.0;
 if (quda_inv_param.dslash_type == QUDA_CLOVER_WILSON_DSLASH){
    if (!QudaInfo::clover_on_device){
       rolex.reset(); rolex.start();
       loadCloverQuda(NULL, NULL, &quda_inv_param);  // allocates space for and compute clover field
       printLaph("...clover term set up and loaded on device");
       printLaph("...clover term should be positive-definite, then clover inverse calculated using Cholesky");
       rolex.stop();
       clovertime=rolex.getTimeInSeconds();
       printLaph(make_str(" Time to set up clover term was ",clovertime," seconds"));
       QudaInfo::clover_on_device=true;}}

     // load the LapH eigenvectors onto host, store pointers
 rolex.reset(); rolex.start();
 int nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 if (!qSmearHandler->queryLaphEigenvectors()){
    errorLaph("Could not load the LapH eigenvectors--computation cannot continue");}
 const vector<LattField>& laphevs(qSmearHandler->getLaphEigenvectors());
 vector<void*> evList(nEigs);
 for (int n=0;n<nEigs;n++){
    evList[n] = (void*)(laphevs[n].getDataConstPtr());}
 rolex.stop();
 double evreadtime=rolex.getTimeInSeconds();
 printLaph(make_str("All Laph eigvecs read in ",evreadtime," seconds\n"));
 rolex.reset(); rolex.start();
 vector<LattField> scaled_evs(laphevs);
 vector<void*> evList_scaled(nEigs);
 for (int n=0;n<nEigs;n++){
    if (insHandler) { 
    	PerambulatorHandler::scale_by_ins_times_coord(scaled_evs[n],insHandler->getInsertion().getTwistVector());
    } else if (sginsHandler) {
	PerambulatorHandler::scale_by_ins_times_coord(scaled_evs[n],sginsHandler->getInsertion().getTwistVector());}
    evList_scaled[n]=(void*)(scaled_evs[n].getDataConstPtr());}
 rolex.stop();
 double evscaledreadtime=rolex.getTimeInSeconds();
 printLaph(make_str("All rescaled Laph eigvecs read in ",evscaledreadtime," seconds\n"));

    // set up any preconditioning in the inverter
 double precondtime=0.0;
 rolex.reset(); rolex.start();
 bool precond=setUpPreconditioning(quda_inv_param);
 rolex.stop();
 precondtime=rolex.getTimeInSeconds();
 if (precond){
    printLaph(make_str(" Time to set up inverter preconditioner was ",precondtime," seconds"));}
 
    // set up Insertion QudaInversionParam
    if( insHandler ) { (insHandler->getInsertion()).setQudaInvertParam(quda_ins_param,*qactionPtr); }

 double srctime=0.0;
 double invtime=0.0;
 double suminstime=0.0;
 double evprojtime=0.0;
 double writetime=0.0;
 int count=0;
 int ncomp=perambComps.computations.size();

 for (list<PerambComputation>::const_iterator it=perambComps.computations.begin();
      it!=perambComps.computations.end();count++,it++){
    printLaph("\n\n *************************************************************");
    printLaph(make_strf(" *\n *  Now starting computation %d (with %d as last):",
               count,ncomp-1));
    printLaph(" *");
    printLaph(" *************************************************************");

    computePerambulators(it->src_time,it->tmin,it->tmax,it->src_lapheigvec_indices,evList,evList_scaled,print_coeffs,extra_soln_check,
                         report_gflops,srctime,invtime,suminstime,evprojtime,writetime,perambComps.useMultiSrcInverter);
    }

 totaltime.stop();
 printLaph("\n\n");
 printLaph("computePerambulators: ran successfully");
 printLaph(make_str("                       Total time = ",totaltime.getTimeInSeconds()," seconds"));
 printLaph(make_str(" Time to load gauge configuration = ",grtime," seconds"));
 printLaph(make_str("            LapH eigvec read time = ",evreadtime," seconds"));
 if (quda_inv_param.dslash_type == QUDA_CLOVER_WILSON_DSLASH){ 
    printLaph(make_str("       Time to set up clover term = ",clovertime," seconds"));}
 printLaph(make_str("         Total source set up time = ",srctime," seconds"));
 if (precond){
    printLaph(make_str("     Inverter preconditioner time = ",precondtime," seconds"));}
 printLaph(make_str("           Total time in inverter = ",invtime," seconds"));
 printLaph(make_str("    Total time applying Insertion = ",suminstime," seconds"));
 printLaph(make_str("Projection onto LapH eigvecs time = ",evprojtime," seconds"));
 printLaph(make_str("          Total file writing time = ",writetime," seconds\n\n"));
}


void PerambulatorHandler::computePerambulators(int src_time, int tmin, int tmax, const set<int>& src_evindices,
                                               const std::vector<void*>& evList, const std::vector<void*>& evList_scaled,
                                               bool print_coeffs, bool extra_soln_check, bool report_gflops, 
                                               double& makesrc_time, double& inv_time, double& ins_time,
                                               double& evproj_time, double& write_time,
                                               bool use_multisrc_inverter)
{
 if (use_multisrc_inverter){
    computePerambulatorsMS(src_time,tmin,tmax,src_evindices,evList,evList_scaled,print_coeffs,extra_soln_check,report_gflops, 
                           makesrc_time,inv_time,ins_time,evproj_time,write_time); }
 else{
    computePerambulatorsSS(src_time,tmin,tmax,src_evindices,evList,print_coeffs,extra_soln_check,report_gflops, 
                           makesrc_time,inv_time,ins_time,evproj_time,write_time); }
}


    //  do inversions for one source time, one requested selection
    //  of source laph-eigenvector dilution indices (all source spins) 
    // VERSION which uses multi-src (MS) inverter

void PerambulatorHandler::computePerambulatorsMS(int src_time, int tmin, int tmax, const set<int>& src_evindices,
                                                 const std::vector<void*>& evList, const std::vector<void*>& evList_scaled,
                                                 bool print_coeffs, bool extra_soln_check, bool report_gflops, 
                                                 double& makesrc_time, double& inv_time, double& sumins_time,
                                                 double& evproj_time, double& write_time)
{
 StopWatch bulova; bulova.start();
 printLaph("\nQuark perambulator computation for one source time,");
 printLaph(" one set of source eigvec indices beginning");
 printLaph(make_strf(" Source time = %d",src_time));

 int Textent = uPtr->getTimeExtent();
 int Lx = LayoutInfo::getLattExtents()[0];
 int Ly = LayoutInfo::getLattExtents()[1];
 int Lz = LayoutInfo::getLattExtents()[2];
 int nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 int minTime=0;
 int maxTime=Textent-1;
 uint nSinkLaphBatch=perambComps.nSinkLaphBatch;
 uint nSinkQudaBatch=perambComps.nSinkQudaBatch;
 uint nEigQudaBatch=perambComps.nEigQudaBatch;
 double soln_rescale=qactionPtr->getSolutionRescaleFactor();


        // allocate space for batched solutions and sources, and make pointers suitable for quda
 vector<int> sinkBatchInds(nSinkLaphBatch);
 vector<int> sinkBatchDoneInds(nSinkLaphBatch);
 vector<int> sinkBatchIndsMap(nSinkLaphBatch);
 vector<LattField> sinkBatchData(nSinkLaphBatch,FieldSiteType::ColorSpinVector);
 vector<LattField> srcBatchData(nSinkLaphBatch,FieldSiteType::ColorSpinVector);
 vector<void*> sinkList(nSinkLaphBatch);
 vector<void*> srcList(nSinkLaphBatch);
 vector<void*> doneList(nSinkLaphBatch);
 for (int iSink=0; iSink<int(nSinkLaphBatch); iSink++){
    sinkList[iSink] = (void*)(sinkBatchData[iSink].getDataPtr());
    srcList[iSink] = (void*)(srcBatchData[iSink].getDataPtr());}
 void** sinks_ptr=(void**)sinkList.data();
 void** srcs_ptr=(void**)srcList.data();
 void** evs_ptr=(void**)evList.data();
 void** evs_scaled_ptr=(void**)evList_scaled.data();
 void** done_sinks_ptr=(void**)doneList.data();

 bulova.stop();
 double srctime=bulova.getTimeInSeconds();
 double invtime=0.0;
 double suminstime=0.0;
 double writetime=0.0;
 double evprojtime=0.0;
 int ninv=src_evindices.size();

     // loop over source laph eigvec indices and source spin

 for (int srcspin=1;srcspin<=int(Nspin);++srcspin){

     // get the file key and open file for writing
  FileKey fkey(src_time,srcspin);
  if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD) == 1) 
   DHputPtr->open(fkey);
  if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID) == 1) 
   DHputPtrSparseGrid->open(fkey);
  if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION) == 1) 
   DHputPtrSummedInsertion->open(fkey);
  if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION) == 1) 
   DHputPtrLeftInsertion->open(fkey);
  if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION) == 1) 
   DHputPtrRightInsertion->open(fkey);
  if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION) == 1) 
   DHputPtrSparseGridRightInsertion->open(fkey);
  if ( (*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION) == 1) 
   DHputPtrSparseGridSummedInsertion->open(fkey);

  int invcount=0;
  int iSinkBatch = 0;

  for (set<int>::const_iterator vt=src_evindices.begin();vt!=src_evindices.end();++vt,++invcount){

    int srcev_ind=*vt;
    printLaph(make_strf("\nStarting source number %d (with %d as last) with source spin %d",invcount,ninv-1,srcspin));
    bool doneflag=true;
    for (int sink_time=0;sink_time<Textent;++sink_time){
    for (int sinkspin=1;sinkspin<=int(Nspin);++sinkspin){
       if (!DHputPtr->queryData(RecordKey(sinkspin,sink_time,srcev_ind))){
          doneflag=false; sinkspin=Nspin+1; sink_time=Textent;}}}

    if ((!fPtr->isModeOverwrite())&&(doneflag)){
       printLaph("warning: these quark sinks already computed...");
       printLaph("  skip re-computing since fileMode not overwrite");}
    else{
        //  initialize source (include gamma_4) in Dirac-Pauli basis, allocate sink
       bulova.reset();bulova.start();
       if  ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)==1){
            make_source(srcBatchData[iSinkBatch],evList_scaled[srcev_ind],src_time,srcspin);
       } else if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION)==1){
 	         make_source(srcBatchData[iSinkBatch],evList_scaled[srcev_ind],src_time,srcspin);
       } else {
            make_source(srcBatchData[iSinkBatch],evList[srcev_ind],src_time,srcspin);
       }
       bulova.stop();
       double addsrctime=bulova.getTimeInSeconds();
       printLaph(make_str(" fermion source set up in Dirac-Pauli basis in ",addsrctime," seconds"));
       srctime+=addsrctime;
       sinkBatchInds[iSinkBatch] = srcev_ind;
       iSinkBatch++; }

       // do the inversions if all sources made
       //   (source and sink should be on host in Dirac-Pauli basis)
    if ((iSinkBatch == int(nSinkLaphBatch))||((invcount==(ninv-1))&&(iSinkBatch>0))) {
       bulova.reset(); bulova.start();
       int nSinks=iSinkBatch;
       iSinkBatch=0;
       quda_inv_param.num_src=nSinks;

       invertMultiSrcQuda(sinks_ptr, srcs_ptr, &quda_inv_param);
 
       bulova.stop(); 
       double addinvtime=bulova.getTimeInSeconds();
       invtime+=addinvtime;
       printLaph(make_strf("Inversions done:  number of iterations = %d",quda_inv_param.iter));
       printLaph(make_str(" Time for ",nSinks," inversions was ",addinvtime," seconds"));
       if (report_gflops) printLaph(make_strf("Multi-source inversion speed was %g gflops",quda_inv_param.gflops));
       bool somefail=false;

       if (!((quda_inv_param.iter>0)&&(quda_inv_param.iter<int(invertPtr->getMaxIterations())))){
          printLaph("\n\nSome inversions in this batch FAILED to converge before max iteration reached");
          somefail=true;
          printLaph("Some solution NOT WRITTEN to file\n\n");}

       int iSinkDone=0;
       for (int iSink=0;iSink<nSinks;++iSink){
          printLaph(make_strf(" Residual for src %d = %g",iSink,quda_inv_param.true_res[iSink]));
          if ((somefail)&&(quda_inv_param.true_res[iSink]>quda_inv_param.tol)){
             printLaph("    This inversion failed to reach required tolerance");}
          else{

             sinkBatchDoneInds[iSinkDone]=sinkBatchInds[iSink];
	          sinkBatchIndsMap[iSinkDone]=iSink;
             doneList[iSinkDone]=(void*)(sinkBatchData[iSink].getDataPtr());
             iSinkDone++;

             // additional checks of the solution

             if (extra_soln_check){
                printLaph("Performing additional check on solution using quda::MatQuda");
                LattField sol_check(FieldSiteType::ColorSpinVector);
                void *sol_checkptr=sol_check.getDataPtr();
                MatQuda(sol_checkptr,sinkList[iSink],&quda_inv_param);
                LaphEnv::compare_latt_fields(sol_check,srcBatchData[iSink]);}
             }}
       nSinks=iSinkDone;


       // carry out projections
       if (nSinks>0){

	       bulova.reset(); bulova.start();
	       int nSinksBatch=std::min(nSinks,int(nSinkQudaBatch));

	       Array<dcmplx> qudaRes(Nspin, Textent, nEigs, nSinks);   // quda laph reversing major order
	       __complex__ double* qudaResPtr = (__complex__ double*)(&qudaRes(0,0,0,0));

	       // sparse grid output to file
			 if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID)==1) {
                gridProject (doneList, sinkBatchDoneInds, nSinks, minTime, maxTime, Textent, print_coeffs, DHputPtrSparseGrid);
					 bulova.stop();
					 double otime=bulova.getTimeInSeconds();
					 printLaph(make_str(" Output of this batch to sparse grid file took ",otime," seconds"));
					 writetime+=otime;
					 bulova.reset(); bulova.start();
				 }
          // sparse grid right insertion output to file
          if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION)==1) {
                gridProject (doneList, sinkBatchDoneInds, nSinks, minTime, maxTime, Textent, print_coeffs, DHputPtrSparseGridRightInsertion);
					 bulova.stop();
					 double otime=bulova.getTimeInSeconds();
					 printLaph(make_str(" Output of this batch to sparse grid file took ",otime," seconds"));
					 writetime+=otime;
					 bulova.reset(); bulova.start();
				 }
          // left insertion output to file
          if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION)==1) {
	          printLaph("projecting batch of solutions onto insertion rescaled LapH eigenvectors");
	          
             laphSinkProject(qudaResPtr, done_sinks_ptr, nSinks, nSinksBatch, evs_scaled_ptr, nEigs, 
			          nEigQudaBatch, &quda_inv_param, LayoutInfo::getRankLattExtents().data());

	          bulova.stop(); 
	          double addevprojtime=bulova.getTimeInSeconds();
	          printLaph(make_str(" this batch of projections onto Laph evs took ",addevprojtime," seconds"));
	          evprojtime+=addevprojtime;

	          // rearrange data then output to file
	          bulova.reset(); bulova.start();
	          for (int iSink=0; iSink<nSinks; ++iSink) {
		          for (int tind=minTime;tind<=maxTime;tind++){
			   			 int t = (tind+Textent)%Textent; 
			          for (int iSpin=0; iSpin<int(Nspin); ++iSpin) {
			   	       vector<dcmplx> quark_sink(nEigs);
			   	       for (int iEv=0; iEv<nEigs; ++iEv) {
			   		       quark_sink[iEv] = soln_rescale*qudaRes(iSpin, t, iEv, iSink);}
		   		       DHputPtrLeftInsertion->putData(RecordKey(iSpin+1,t,sinkBatchDoneInds[iSink]),quark_sink);
	   			       if (print_coeffs){
		   			       printLaph(make_strf("srcev_index = %d, spin = %d, time = %d",sinkBatchDoneInds[iSink],iSpin+1,t));
			   		       for (int n=0;n<nEigs;n++){
				   		       printLaph(make_strf("coef for eigenlevel %d = (%14.8f, %14.8f)",
					   				       n,real(quark_sink[n]),imag(quark_sink[n])));}}}}}
					          bulova.stop();
					          double otime=bulova.getTimeInSeconds();
					          printLaph(make_str(" Output of this batch to file took ",otime," seconds"));
					          writetime+=otime;
             }
	       // do the Laph - projections and standard output
          if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD)==1){
	          printLaph("projecting batch of solutions onto LapH eigenvectors");
	          laphSinkProject(qudaResPtr, done_sinks_ptr, nSinks, nSinksBatch, evs_ptr, nEigs, 
			          nEigQudaBatch, &quda_inv_param, LayoutInfo::getRankLattExtents().data());

	          bulova.stop(); 
	          double addevprojtime=bulova.getTimeInSeconds();
	          printLaph(make_str(" this batch of projections onto Laph evs took ",addevprojtime," seconds"));
	          evprojtime+=addevprojtime;

	          // rearrange data then output to file
	          bulova.reset(); bulova.start();
	          for (int iSink=0; iSink<nSinks; ++iSink) {
		          for (int tind=minTime;tind<=maxTime;tind++){
			   			 int t = (tind+Textent)%Textent; 
			          for (int iSpin=0; iSpin<int(Nspin); ++iSpin) {
			   	       vector<dcmplx> quark_sink(nEigs);
			   	       for (int iEv=0; iEv<nEigs; ++iEv) {
			   		       quark_sink[iEv] = soln_rescale*qudaRes(iSpin, t, iEv, iSink);}
		   		       DHputPtr->putData(RecordKey(iSpin+1,t,sinkBatchDoneInds[iSink]),quark_sink);
	   			       if (print_coeffs){
		   			       printLaph(make_strf("srcev_index = %d, spin = %d, time = %d",sinkBatchDoneInds[iSink],iSpin+1,t));
			   		       for (int n=0;n<nEigs;n++){
				   		       printLaph(make_strf("coef for eigenlevel %d = (%14.8f, %14.8f)",
					   				       n,real(quark_sink[n]),imag(quark_sink[n])));}}}}}
					          bulova.stop();
					          double otime=bulova.getTimeInSeconds();
					          printLaph(make_str(" Output of this batch to file took ",otime," seconds"));
					          writetime+=otime;
             }
          // right insertion output to file
          if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)==1){
	          printLaph("projecting batch of solutions onto LapH eigenvectors");
	          laphSinkProject(qudaResPtr, done_sinks_ptr, nSinks, nSinksBatch, evs_ptr, nEigs, 
			          nEigQudaBatch, &quda_inv_param, LayoutInfo::getRankLattExtents().data());

	          bulova.stop(); 
	          double addevprojtime=bulova.getTimeInSeconds();
	          printLaph(make_str(" this batch of projections onto Laph evs took ",addevprojtime," seconds"));
	          evprojtime+=addevprojtime;

	          // rearrange data then output to file
	          bulova.reset(); bulova.start();
	          for (int iSink=0; iSink<nSinks; ++iSink) {
		          for (int tind=minTime;tind<=maxTime;tind++){
			   			 int t = (tind+Textent)%Textent; 
			          for (int iSpin=0; iSpin<int(Nspin); ++iSpin) {
			   	       vector<dcmplx> quark_sink(nEigs);
			   	       for (int iEv=0; iEv<nEigs; ++iEv) {
			   		       quark_sink[iEv] = soln_rescale*qudaRes(iSpin, t, iEv, iSink);}
		   		       DHputPtrRightInsertion->putData(RecordKey(iSpin+1,t,sinkBatchDoneInds[iSink]),quark_sink);
	   			       if (print_coeffs){
		   			       printLaph(make_strf("srcev_index = %d, spin = %d, time = %d",sinkBatchDoneInds[iSink],iSpin+1,t));
			   		       for (int n=0;n<nEigs;n++){
				   		       printLaph(make_strf("coef for eigenlevel %d = (%14.8f, %14.8f)",
					   				       n,real(quark_sink[n]),imag(quark_sink[n])));}}}}}
					          bulova.stop();
					          double otime=bulova.getTimeInSeconds();
					          printLaph(make_str(" Output of this batch to file took ",otime," seconds"));
					          writetime+=otime;
             }

             // carry out summed insertion 
          if (((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)==1)||
              ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION)==1)) {
             bulova.reset(); bulova.start();
             for (int iSink=0; iSink<nSinks; ++iSink){
                // applyInsertion(srcBatchData[iSink],sinkBatchData[sinkBatchIndsMap[iSink]]); // Insertion-Multiplication is performed on CPU
		          MatQuda(srcList[iSink],doneList[iSink],&quda_ins_param);
               }
             bulova.stop(); 
             double addsuminstime=bulova.getTimeInSeconds();
             suminstime+=addsuminstime;
             printLaph(make_str(" Time for ",nSinks," insertions was ",addsuminstime," seconds"));
             bulova.reset(); bulova.start();
            
             // carry out second inversion 
             quda_inv_param.num_src=nSinks;
             invertMultiSrcQuda(sinks_ptr, srcs_ptr, &quda_inv_param);  
            
             bulova.stop(); 
             double addinvtime=bulova.getTimeInSeconds();
             invtime+=addinvtime;
             printLaph(make_strf("Inversions done:  number of iterations = %d",quda_inv_param.iter));
             printLaph(make_str(" Time for ",nSinks," inversions was ",addinvtime," seconds"));
             if (report_gflops) printLaph(make_strf("Multi-source inversion speed was %g gflops",quda_inv_param.gflops));
             bool somefail=false;

             if (!((quda_inv_param.iter>0)&&(quda_inv_param.iter<int(invertPtr->getMaxIterations())))){
                printLaph("\n\nSome inversions in this batch FAILED to converge before max iteration reached");
                somefail=true; 
                printLaph("Some solution NOT WRITTEN to file\n\n");}

             int iSinkDone=0;
             for (int iSink=0;iSink<nSinks;++iSink){
                printLaph(make_strf(" Residual for src %d = %g",iSink,quda_inv_param.true_res[iSink]));
                if ((somefail)&&(quda_inv_param.true_res[iSink]>quda_inv_param.tol)){
                   printLaph("    This inversion failed to reach required tolerance");}
                else{

                   sinkBatchDoneInds[iSinkDone]=sinkBatchInds[iSink];
                   doneList[iSinkDone]=(void*)(sinkBatchData[iSink].getDataPtr());
                   iSinkDone++;

               // additional checks of the solution

                   if (extra_soln_check){
                      printLaph("Performing additional check on solution using quda::MatQuda");
                      LattField sol_check(FieldSiteType::ColorSpinVector);
                      void *sol_checkptr=sol_check.getDataPtr();
                      MatQuda(sol_checkptr,sinkList[iSink],&quda_inv_param);
                      LaphEnv::compare_latt_fields(sol_check,srcBatchData[iSink]);}
                  }}
             nSinks=iSinkDone;
            
             if (nSinks>0){

                bulova.reset(); bulova.start();
                int nSinksBatch=std::min(nSinks,int(nSinkQudaBatch));

               // sparse grid summed insertion output to file
                if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION)==1) {
                   gridProject (doneList, sinkBatchDoneInds, nSinks, minTime, maxTime, Textent, print_coeffs, DHputPtrSparseGridSummedInsertion);
					    bulova.stop();
					    double otime=bulova.getTimeInSeconds();
					    printLaph(make_str(" Output of this batch to sparse grid insertion file took ",otime," seconds"));
					    writetime+=otime;
					    bulova.reset(); bulova.start();
				       } // end sparse grid summed insertion output to file 

                // do Laph - projections and summed insertion output to file
                if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)==1) {
                   printLaph("projecting batch of insertion solutions onto LapH eigenvectors");
                   laphSinkProject(qudaResPtr, done_sinks_ptr, nSinks, nSinksBatch, evs_ptr, nEigs, 
                         nEigQudaBatch, &quda_inv_param, LayoutInfo::getRankLattExtents().data());
 
                   bulova.stop(); 
                   double addevprojtime=bulova.getTimeInSeconds();
                   printLaph(make_str(" this batch of projections onto Laph evs took ",addevprojtime," seconds"));
                   evprojtime+=addevprojtime;
 
                   // rearrange data then output to file
                   bulova.reset(); bulova.start();
                   for (int iSink=0; iSink<nSinks; ++iSink) {
                      for (int tind=minTime;tind<=maxTime;tind++){
                            int t = (tind+Textent)%Textent; 
                         for (int iSpin=0; iSpin<int(Nspin); ++iSpin) {
                            vector<dcmplx> quark_sink(nEigs);
                            for (int iEv=0; iEv<nEigs; ++iEv) {
                               quark_sink[iEv] = soln_rescale*qudaRes(iSpin, t, iEv, iSink);}
                            DHputPtrSummedInsertion->putData(RecordKey(iSpin+1,t,sinkBatchDoneInds[iSink]),quark_sink);
                            if (print_coeffs){
                               printLaph(make_strf("srcev_index = %d, spin = %d, time = %d",sinkBatchDoneInds[iSink],iSpin+1,t));
                               for (int n=0;n<nEigs;n++){
                                  printLaph(make_strf("coef for eigenlevel %d = (%14.8f, %14.8f)",
                                           n,real(quark_sink[n]),imag(quark_sink[n])));}}}}}
                               bulova.stop();
                               double otime=bulova.getTimeInSeconds();
                               printLaph(make_str(" Output of this batch to file took ",otime," seconds"));
                               writetime+=otime;
                   } // end summed insertion projection output to file
             } // n_Sinks after second inversion > 0
         } // carry out summed insertion
      } // n_Sinks after first inversion > 0
   }  // batch end

    if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD)==1) DHputPtr->flush();
	 if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID)==1) DHputPtrSparseGrid->flush();
    if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)==1) DHputPtrSummedInsertion->flush();
    if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION)==1) DHputPtrLeftInsertion->flush();
    if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)==1) DHputPtrRightInsertion->flush();
    if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION)==1) DHputPtrSparseGridRightInsertion->flush();
    if ((*pmapPtr).count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION)==1) DHputPtrSparseGridSummedInsertion->flush();
	}}   // src_ind, src_spin loop end

 inv_time+=invtime;
 sumins_time+=suminstime;
 makesrc_time+=srctime; 
 evproj_time+=evprojtime; 
 write_time+=writetime;
}

    //  do inversions for one source time, one requested selection
    //  of source laph-eigenvector dilution indices (all source spins) 
    //  VERSION which does not use multi-src inverter, but single src (SS)

void PerambulatorHandler::computePerambulatorsSS(int src_time, int tmin, int tmax, const set<int>& src_evindices,
                                                 const std::vector<void*>& evList, bool print_coeffs,
                                                 bool extra_soln_check, bool report_gflops, 
                                                 double& makesrc_time, double& inv_time, double& ins_time,
                                                 double& evproj_time, double& write_time)
{
 StopWatch bulova; bulova.start();
 printLaph("\nQuark perambulator computation for one source time,");
 printLaph(" one set of source eigvec indices beginning");
 printLaph(make_strf(" Source time = %d, min = %d, max = %d",src_time,tmin,tmax));

 int Textent = uPtr->getTimeExtent();
 int nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 int minTime=tmin;
 int maxTime=tmax;
 uint nSinkLaphBatch=perambComps.nSinkLaphBatch;
 uint nSinkQudaBatch=perambComps.nSinkQudaBatch;
 uint nEigQudaBatch=perambComps.nEigQudaBatch;

 double soln_rescale=qactionPtr->getSolutionRescaleFactor();

 // allocate space for batched solutions and make pointers suitable for quda
 int iSinkBatch = 0;
 vector<int> sinkBatchInds(nSinkLaphBatch);
 vector<LattField> sinkBatchData(nSinkLaphBatch,FieldSiteType::ColorSpinVector);
 vector<void*> sinkList(nSinkLaphBatch);
 for (int iSink=0; iSink<int(nSinkLaphBatch); iSink++){
    sinkList[iSink] = (void*)(sinkBatchData[iSink].getDataPtr());}
 void** sinks_ptr=(void**)sinkList.data();
 void** evs_ptr=(void**)evList.data();

 bulova.stop();
 double srctime=bulova.getTimeInSeconds();
 double invtime=0.0;
 double suminstime=0.0;
 double writetime=0.0;
 double evprojtime=0.0;
 int ninv=src_evindices.size();

 if (sgHandler)
	 printLaph("\nIn computeSS: sparse grid enabled.");
 else 
	 printLaph("\nIn computeSS: sparse grid disabled.");
 if (insHandler)
	 printLaph("\nIn computeSS: insertion enabled.");
 else
         printLaph("\nIn computeSS: insertion disabled.");
 if (sginsHandler)
         printLaph("\nIn computeSS: sparse grid projected insertion enabled.");
 else
         printLaph("\nIn computeSS: sparse grid projected insertion disabled.");
   
 // loop over source laph eigvec indices and source spin
 
 for (int srcspin=1;srcspin<=int(Nspin);++srcspin){

     // get the file key and open file for writing
  FileKey fkey(src_time,srcspin);
  DHputPtr->open(fkey);
  if (sgHandler) 
		DHputPtrSparseGrid->open(fkey);
  if (insHandler) DHputPtrSummedInsertion->open(fkey);
  if (sginsHandler) DHputPtrSparseGridSummedInsertion->open(fkey);

  int invcount=0;

  for (set<int>::const_iterator vt=src_evindices.begin();vt!=src_evindices.end();++vt,++invcount){

    int srcev_ind=*vt;
    printLaph(make_strf("\nStarting inversion number %d (with %d as last) with source spin %d",invcount,ninv-1,srcspin));
    bool doneflag=true;
    for (int sink_time=0;sink_time<Textent;++sink_time){
    for (int sinkspin=1;sinkspin<=int(Nspin);++sinkspin){
       if (!DHputPtr->queryData(RecordKey(sinkspin,sink_time,srcev_ind))){
          doneflag=false; sinkspin=Nspin+1; sink_time=Textent;}}}

    if ((!fPtr->isModeOverwrite())&&(doneflag)){
       printLaph("warning: these quark sinks already computed...");
       printLaph("  skip re-computing since fileMode not overwrite");}
    else{

        //  initialize source (include gamma_4) in Dirac-Pauli basis, allocate sink
    bulova.reset();bulova.start();
    LattField ferm_src;
    make_source(ferm_src,evList[srcev_ind],src_time,srcspin);
    bulova.stop();
    double addsrctime=bulova.getTimeInSeconds();
    printLaph(make_str(" fermion source set up in Dirac-Pauli basis in ",addsrctime," seconds"));
    srctime+=addsrctime;

    // now do the inversion (source and sink should be on host in Dirac-Pauli basis)
    bulova.reset(); bulova.start();
    void *spinor_src = (void*)(ferm_src.getDataPtr());
    void *spinor_snk = (void*)(sinkBatchData[iSinkBatch].getDataPtr());

    invertQuda(spinor_snk, spinor_src, &quda_inv_param);

    bulova.stop(); 
    double addinvtime=bulova.getTimeInSeconds();
    invtime+=addinvtime;

    printLaph(make_strf("Inversion done:  number of iterations = %d",quda_inv_param.iter));
    printLaph(make_strf(" Residual = %g",quda_inv_param.true_res[0]));
    printLaph(make_str("  Inversion time was ",addinvtime," seconds"));
    if (report_gflops) printLaph(make_strf("Inversion speed was %g gflops",quda_inv_param.gflops));
    //printLaph(make_str("  Inverter used was ",invertPtr->getName()));

    bool success=((quda_inv_param.iter>0)&&(quda_inv_param.iter<int(invertPtr->getMaxIterations())));
    if (quda_inv_param.true_res[0]>quda_inv_param.tol){ 
       printLaph("Inverter failed to reach required tolerance");}

       // additional checks of the solution

    if (extra_soln_check){
       printLaph("Performing additional check on solution using quda::MatQuda");
       LattField sol_check(FieldSiteType::ColorSpinVector);
       void *sol_checkptr=sol_check.getDataPtr();
       MatQuda(sol_checkptr,spinor_snk,&quda_inv_param);
       LaphEnv::compare_latt_fields(sol_check,ferm_src);} 

    if (success){
       sinkBatchInds[iSinkBatch] = srcev_ind;
       iSinkBatch++;}
    else{
       printLaph("\n\nInversion FAILED to converge before max iteration reached");
       printLaph("Solution NOT WRITTEN to file\n\n");}
    }

       // carry out projections
    if ((iSinkBatch == int(nSinkLaphBatch))||((invcount==(ninv-1))&&(iSinkBatch>0))) {

       bulova.reset(); bulova.start();
       int nSinks=iSinkBatch;
       int nSinksBatch=std::min(nSinks,int(nSinkQudaBatch));

       Array<dcmplx> qudaRes(Nspin, Textent, nEigs, nSinks);   // quda laph reversing major order
       __complex__ double* qudaResPtr = (__complex__ double*)(&qudaRes(0,0,0,0));

			 if (sgHandler) {
				 gridProject (sinkList, sinkBatchInds, nSinks, minTime, maxTime, Textent, print_coeffs, DHputPtrSparseGrid);
				 bulova.stop();
				 double otime=bulova.getTimeInSeconds();
				 printLaph(make_str(" Output of this batch to sparse grid file took ",otime," seconds"));
				 writetime+=otime;
				 bulova.reset(); bulova.start();
			 }
			 // do the projections
			 printLaph("projecting batch of solutions onto LapH eigenvectors");
			 laphSinkProject(qudaResPtr, sinks_ptr, nSinks, nSinksBatch, evs_ptr, nEigs, 
					 nEigQudaBatch, &quda_inv_param, LayoutInfo::getRankLattExtents().data());

			 bulova.stop(); 
			 double addevprojtime=bulova.getTimeInSeconds();
			 printLaph(make_str(" this batch of projections onto Laph evs took ",addevprojtime," seconds"));
			 evprojtime+=addevprojtime;

			 // rearrange data then output to file
			 bulova.reset(); bulova.start();
			 for (int iSink=0; iSink<nSinks; ++iSink) {
				 for (int tind=minTime;tind<=maxTime;tind++){
					 int t = (tind+Textent)%Textent; 
					 for (int iSpin=0; iSpin<int(Nspin); ++iSpin) {  
						 vector<dcmplx> quark_sink(nEigs);
						 for (int iEv=0; iEv<nEigs; ++iEv) {
							 quark_sink[iEv] = soln_rescale*qudaRes(iSpin, t, iEv, iSink);}
						 DHputPtr->putData(RecordKey(iSpin+1,t,sinkBatchInds[iSink]),quark_sink);
						 if (print_coeffs){
							 printLaph(make_strf("srcev_index = %d, spin = %d, time = %d",sinkBatchInds[iSink],iSpin+1,t));
							 for (int n=0;n<nEigs;n++){
								 printLaph(make_strf("coef for eigenlevel %d = (%14.8f, %14.8f)",
											 n,real(quark_sink[n]),imag(quark_sink[n])));}}}}}
			 bulova.stop();
			 double otime=bulova.getTimeInSeconds();
			 printLaph(make_str(" Output of this batch to file took ",otime," seconds"));
			 writetime+=otime;
	   
	  // carry out insertion 
          if (insHandler) {
	     LattField ferm_src(FieldSiteType::ColorSpinVector);
	     void *spinor_src = (void*)(ferm_src.getDataPtr());
	     int iSinkBatchtwo = 0;
	     vector<int> sinkBatchIndstwo(iSinkBatch);
             for (int iSink=0; iSink<iSinkBatch; ++iSink){
                bulova.reset(); bulova.start();
		// applyInsertion(srcBatchData[iSink],sinkBatchData[sinkBatchIndsMap[iSink]]); // Insertion-Multiplication is performed on CPU
                void *spinor_snk = (void*)(sinkBatchData[iSink].getDataPtr()); 
		MatQuda(spinor_src,spinor_snk,&quda_ins_param);
                bulova.stop();
                double addsuminstime=bulova.getTimeInSeconds();
                suminstime+=addsuminstime;
                printLaph(make_str(" Insertion Time was ",addsuminstime," seconds"));
                bulova.reset(); bulova.start();
		
		// carry out second inversion
	        invertQuda(spinor_snk, spinor_src, &quda_inv_param);

	        bulova.stop();
      	        double addinvtime=bulova.getTimeInSeconds();
		invtime+=addinvtime;

		printLaph(make_strf("Inversion done:  number of iterations = %d",quda_inv_param.iter));
	        printLaph(make_strf(" Residual = %g",quda_inv_param.true_res[0]));
     	        printLaph(make_str("  Inversion time was ",addinvtime," seconds"));
		if (report_gflops) printLaph(make_strf("Inversion speed was %g gflops",quda_inv_param.gflops));
		//printLaph(make_str("  Inverter used was ",invertPtr->getName()));

		bool success=((quda_inv_param.iter>0)&&(quda_inv_param.iter<int(invertPtr->getMaxIterations())));
		if (quda_inv_param.true_res[0]>quda_inv_param.tol){
		    printLaph("Inverter failed to reach required tolerance");}

		    // additional checks of the solution

		if (extra_soln_check){
		   printLaph("Performing additional check on solution using quda::MatQuda");
		   LattField sol_check(FieldSiteType::ColorSpinVector);
		   void *sol_checkptr=sol_check.getDataPtr();
		   MatQuda(sol_checkptr,spinor_snk,&quda_inv_param);
		   LaphEnv::compare_latt_fields(sol_check,ferm_src);}

		if (success){
		   sinkBatchIndstwo[iSinkBatchtwo] = sinkBatchInds[iSink];
		   iSinkBatchtwo++;
		  }
		else{
		   printLaph("\n\nInversion FAILED to converge before max iteration reached");
		   printLaph("Solution NOT WRITTEN to file\n\n");}
	     	}

	        // carry out projections
                         if (sginsHandler) {
                                 gridProject (sinkList, sinkBatchIndstwo, iSinkBatchtwo, minTime, maxTime, Textent, print_coeffs, DHputPtrSparseGridSummedInsertion);
                                 bulova.stop();
                                 double otime=bulova.getTimeInSeconds();
                                 printLaph(make_str(" Output of this batch to sparse grid insertion file took ",otime," seconds"));
                                 writetime+=otime;
                                 bulova.reset(); bulova.start();
                         }
                         // do the projections
                         printLaph("projecting batch of insertion solutions onto LapH eigenvectors");
                         laphSinkProject(qudaResPtr, sinks_ptr, nSinks, nSinksBatch, evs_ptr, nEigs,
                                         nEigQudaBatch, &quda_inv_param, LayoutInfo::getRankLattExtents().data());

                         bulova.stop();
                         double addevprojtime=bulova.getTimeInSeconds();
                         printLaph(make_str(" this batch of projections onto Laph evs took ",addevprojtime," seconds"));
                         evprojtime+=addevprojtime;

                         // rearrange data then output to file
                         bulova.reset(); bulova.start();
                         for (int iSink=0; iSink<nSinks; ++iSink) {
                                 for (int tind=minTime;tind<=maxTime;tind++){
                                         int t = (tind+Textent)%Textent;
                                         for (int iSpin=0; iSpin<int(Nspin); ++iSpin) {
                                                 vector<dcmplx> quark_sink(nEigs);
                                                 for (int iEv=0; iEv<nEigs; ++iEv) {
                                                         quark_sink[iEv] = soln_rescale*qudaRes(iSpin, t, iEv, iSink);}
                                                 DHputPtrSummedInsertion->putData(RecordKey(iSpin+1,t,sinkBatchInds[iSink]),quark_sink);
                                                 if (print_coeffs){
                                                         printLaph(make_strf("srcev_index = %d, spin = %d, time = %d",sinkBatchInds[iSink],iSpin+1,t));
                                                         for (int n=0;n<nEigs;n++){
                                                                 printLaph(make_strf("coef for eigenlevel %d = (%14.8f, %14.8f)",
                                                                                         n,real(quark_sink[n]),imag(quark_sink[n])));}}}}}
                         bulova.stop();
                         double otime=bulova.getTimeInSeconds();
                         printLaph(make_str(" Output of this batch to file took ",otime," seconds\n\n"));
                         writetime+=otime;
    			 }

			 iSinkBatch = 0;}    // batch end
		DHputPtr->flush();
	  if (sgHandler)	
			DHputPtrSparseGrid->flush();
	  if (insHandler)
         		DHputPtrSummedInsertion->flush();
     	  if (sginsHandler)
         		DHputPtrSparseGridSummedInsertion->flush();
    	}}   // src_ind, src_spin loop end

 inv_time+=invtime;
 ins_time+=suminstime;
 makesrc_time+=srctime; 
 evproj_time+=evprojtime; 
 write_time+=writetime;
}


// ***************************************************************************

// Creates the source (multiplied by gamma_4) in Dirac-Pauli basis
//   src_spin here is 1,2,3,4

void PerambulatorHandler::make_source(LattField& ferm_src, const void* ev_src_ptr, 
		int src_time, int src_spin)
{
	printLaph(" Making source for this inversion...");
	ferm_src.reset(FieldSiteType::ColorSpinVector);
	bool dp=(ferm_src.bytesPerWord()==sizeof(std::complex<double>));
	// initialize source field to zero
	int loc_nsites=LayoutInfo::getRankLatticeNumSites();
	int ncmplx_per_site=ferm_src.elemsPerSite();
	int ncmplx=ncmplx_per_site*loc_nsites;
	int cbytes;
	dcmplx zrhodp;
	fcmplx zrhosp;
	char *zrho;
	if (dp){
		double* z0=reinterpret_cast<double*>(ferm_src.getDataPtr());
		std::fill(z0,z0+2*ncmplx,0.0);
		cbytes=sizeof(std::complex<double>);
		zrho=reinterpret_cast<char*>(&zrhodp);}
	else{
		float* z0=reinterpret_cast<float*>(ferm_src.getDataPtr());
		std::fill(z0,z0+2*ncmplx,0.0);
		cbytes=sizeof(std::complex<float>);
		zrho=reinterpret_cast<char*>(&zrhosp);}

	int loc_npsites=loc_nsites/2;
	int start_parity=LayoutInfo::getMyStartParity();
	int mytmin=LayoutInfo::getMyCommCoords()[3]*LayoutInfo::getRankLattExtents()[3];
	int mytmax=mytmin+LayoutInfo::getRankLattExtents()[3]-1;
	int tstride=LayoutInfo::getRankLattExtents()[0]*LayoutInfo::getRankLattExtents()[1]
		*LayoutInfo::getRankLattExtents()[2];
	int incx=FieldNcolor;
	int incy=FieldNcolor*FieldNspin;

	if ((src_time>=mytmin)&&(src_time<=mytmax)){
		int tloc=src_time-mytmin;
		int parshift=loc_npsites*((start_parity+tloc)%2);
		int start1=((tstride*tloc)/2) + parshift;
		int stop1=((1+tstride*(tloc+1))/2) + parshift;
		int n1=stop1-start1;
		parshift=loc_npsites*((start_parity+1+tloc)%2);
		int start2=((1+tstride*tloc)/2) + parshift;
		int stop2=((tstride*(tloc+1))/2) + parshift;
		int n2=stop2-start2;
		int xstart1=start1*incx*cbytes;
		int xstart2=start2*incx*cbytes;
		char* ystart1=ferm_src.getDataPtr()+start1*incy*cbytes;
		char* ystart2=ferm_src.getDataPtr()+start2*incy*cbytes;

		const char* x0=reinterpret_cast<const char*>(ev_src_ptr);
		zrhodp=dcmplx(1.0,0.0);
		if (src_spin>2){ zrhodp=-zrhodp;}  // multiply by gamma_4
    if (!dp){
       zrhosp=complex<float>(real(zrhodp),imag(zrhodp));}
    const char* x1=x0+xstart1;
    const char* x2=x0+xstart2;
    char* y1=ystart1+(src_spin-1)*incx*cbytes;
    char* y2=ystart2+(src_spin-1)*incx*cbytes;
    for (int c=0;c<FieldNcolor;++c){
       if (dp){
          cblas_zaxpy(n1,(dcmplx*)(zrho),(dcmplx*)(x1),incx,(dcmplx*)(y1),incy);
          cblas_zaxpy(n2,(dcmplx*)(zrho),(dcmplx*)(x2),incx,(dcmplx*)(y2),incy);}
       else{
          cblas_caxpy(n1,(fcmplx*)(zrho),(fcmplx*)(x1),incx,(fcmplx*)(y1),incy);
          cblas_caxpy(n2,(fcmplx*)(zrho),(fcmplx*)(x2),incx,(fcmplx*)(y2),incy);}
       x1+=cbytes; y1+=cbytes; x2+=cbytes; y2+=cbytes;}}
 printLaph("Source for this inversion created");
}

template <typename T>   // T = double or float
void PerambulatorHandler::scale_by_ins_times_coord(T* data, int ncmplx_per_site,const std::vector<int>& coeff)
{
 const std::vector<int>& ext=LayoutInfo::getRankLattExtents();
 const std::vector<int>& cc=LayoutInfo::getMyCommCoords();
 const int nsites=LayoutInfo::getRankLatticeNumSites();
 const int sp=LayoutInfo::getMyStartParity();
 const int X0=cc[0]*ext[0], Y0=cc[1]*ext[1], Z0=cc[2]*ext[2]; // global coordinates of rank's origin
 
 // loop over all local lattice coordinates
 std::vector<int> loc(4);
 for (loc[3]=0;loc[3]<ext[3];++loc[3])
 for (loc[2]=0;loc[2]<ext[2];++loc[2])
 for (loc[1]=0;loc[1]<ext[1];++loc[1])
 for (loc[0]=0;loc[0]<ext[0];++loc[0]){
    long idx = (long)((loc[0]+ext[0]*(loc[1]+ext[1]*(loc[2]+ext[2]*loc[3])))/2)
          + (long)(nsites/2)*((sp+loc[0]+loc[1]+loc[2]+loc[3])%2); 
    const T pre=T(coeff[0]*(X0+loc[0])+coeff[1]*(Y0+loc[1])+coeff[2]*(Z0+loc[2])); // twist * global coordinates
    T* p=data+2*idx*ncmplx_per_site;
    for (int k=0;k<ncmplx_per_site;++k){
       T re=p[2*k], im=p[2*k+1];
       p[2*k]  =-pre*im;
       p[2*k+1]= pre*re;}}
}

void PerambulatorHandler::scale_by_ins_times_coord(LattField& field,const std::vector<int>& coeff)
{
 int n=field.elemsPerSite();
 bool dp=(field.bytesPerWord()==sizeof(std::complex<double>));
 if (dp) scale_by_ins_times_coord(reinterpret_cast<double*>(field.getDataPtr()),n,coeff);
 else    scale_by_ins_times_coord(reinterpret_cast<float*>(field.getDataPtr()),n,coeff);
}

template <typename H>
void PerambulatorHandler::gridProject (const vector<void*>& _doneList, const vector<int>&  _sinkBatchDoneInds, int _nSinks, int _minTime, int _maxTime, int _Textent, 
		bool _print_coeffs, DataPutHandlerMF<H,PerambulatorHandler::FileKey,PerambulatorHandler::RecordKey,PerambulatorHandler::DataType>* _DHputPtr){
 int nColor = 3; 
 int nGridPoints = sgHandler->getGrid().getNGridPoints(); // not necessary to implement sginsHandler version, since the are based on the same RandomSparseGridInfo
 vector<vector<int>> offsets; sgHandler->getGrid().generateOffsets(offsets); 
 printLaph(" Outputting to sparse grid file"); 
 for (int iSink=0; iSink<_nSinks; ++iSink) {
	 const char* field_start=reinterpret_cast<char*>(_doneList[iSink]);
	 size_t site_bytes=2*sizeof(double)*12;
	 size_t word_bytes=2*sizeof(double);
	 if (word_bytes!=sizeof(complex<double>))
		 throw logic_error("only implemented for double precision"); 
	 size_t colvec_bytes=word_bytes*FieldNcolor;  
	 for (int tind=_minTime;tind<=_maxTime;tind++){
		 int t = (tind+_Textent)%_Textent; 
		 const auto& local_offsets =  
				 (*sgHandler).getGrid().getLocalGridPoints(offsets[t],t);
		 vector<dcmplx> all_spin_quark_sink(nColor*nGridPoints*Nspin,0.0);
		 for (const auto& offset : local_offsets) { 
			 size_t local_offset=site_bytes*offset.local_offset;
			 for (int iSpin=0; iSpin<int(Nspin); ++iSpin) {
				 const char* get_ptr = field_start+local_offset+iSpin*colvec_bytes;
				 size_t global_offset=colvec_bytes*(offset.global_offset + iSpin*nGridPoints);
				 char* dest_ptr = reinterpret_cast<char*>(all_spin_quark_sink.data())+global_offset; 
				 memcpy(dest_ptr,get_ptr,colvec_bytes); 
			 }
		 }
		 MPI_Allreduce(MPI_IN_PLACE, all_spin_quark_sink.data(), 2*all_spin_quark_sink.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		 for (int iSpin=0; iSpin<int(Nspin); ++iSpin) {
		    vector<dcmplx> quark_sink(all_spin_quark_sink.cbegin()+iSpin*nColor*nGridPoints,all_spin_quark_sink.cbegin()+(iSpin+1)*nColor*nGridPoints);
		    _DHputPtr->putData(RecordKey(iSpin+1,t,_sinkBatchDoneInds[iSink]),quark_sink);
		    if (_print_coeffs){
		   	 printLaph(make_strf("srcev_index = %d, spin = %d, time = %d",_sinkBatchDoneInds[iSink],iSpin+1,t));
		   	 for (int n=0;n<(nColor*nGridPoints);n++){
		   		 printLaph(make_strf("component for color/space component %d = (%14.8f, %14.8f)",n,real(quark_sink[n]),imag(quark_sink[n])));
		   	 }
		    }
		 }
	 }
 }
}

const PerambulatorHandler::DataType& PerambulatorHandler::getData(
                             int snk_time, int snk_spin, int src_time, int src_spin,
                             int src_eigvec_index) const
{
 if (!isInfoSet()){
    errorLaph("Cannot getData in PerambulatorHandler unless info is set");}
 if (mode==Compute){
    errorLaph("Cannot getData in PerambulatorHandler unless not in compute mode");}
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(src_eigvec_index<0)||(src_eigvec_index>=nEigs)){
    errorLaph("Invalid parameters passed to getData in PerambulatorHandler");}
 FileKey fkey(src_time,src_spin);
 RecordKey rkey(snk_spin,snk_time,src_eigvec_index);
 return DHgetPtr->getData(fkey,rkey);
}

const PerambulatorHandler::DataType& PerambulatorHandler::getDataSparseGrid(
                             int snk_time, int snk_spin, int src_time, int src_spin,
                             int src_eigvec_index) const
{
 if (!sgHandler){
    errorLaph("Cannot getDataSparseGrid in PerambulatorHandler: no sparse grid specified");}
 if (!isInfoSet()){
    errorLaph("Cannot getData in PerambulatorHandler unless info is set");}
 if (mode==Compute){
    errorLaph("Cannot getData in PerambulatorHandler unless not in compute mode");}
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(src_eigvec_index<0)||(src_eigvec_index>=nEigs)){
    errorLaph("Invalid parameters passed to getData in PerambulatorHandler");}
 FileKey fkey(src_time,src_spin);
 RecordKey rkey(snk_spin,snk_time,src_eigvec_index);
 return DHgetPtrSparseGrid->getData(fkey,rkey);
}

const PerambulatorHandler::DataType& PerambulatorHandler::getDataSummedInsertion(
                             int snk_time, int snk_spin, int src_time, int src_spin,
                             int src_eigvec_index) const
{
 if (!insHandler){
    errorLaph("Cannot getDataInsertion in PerambulatorHandler: no insertion specified");}
 if (!isInfoSet()){
    errorLaph("Cannot getData in PerambulatorHandler unless info is set");}
 if (mode==Compute){
    errorLaph("Cannot getData in PerambulatorHandler unless not in compute mode");}
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(src_eigvec_index<0)||(src_eigvec_index>=nEigs)){
    errorLaph("Invalid parameters passed to getData in PerambulatorHandler");}
 FileKey fkey(src_time,src_spin);
 RecordKey rkey(snk_spin,snk_time,src_eigvec_index);
 return DHgetPtrSummedInsertion->getData(fkey,rkey);
}

const PerambulatorHandler::DataType& PerambulatorHandler::getDataSparseGridSummedInsertion(
                             int snk_time, int snk_spin, int src_time, int src_spin,
                             int src_eigvec_index) const
{
 if (!sginsHandler){
    errorLaph("Cannot getDataSparseGridSummedInsertion in PerambulatorHandler: no sparse grid projected insertion specified");}
 if (!isInfoSet()){
    errorLaph("Cannot getData in PerambulatorHandler unless info is set");}
 if (mode==Compute){
    errorLaph("Cannot getData in PerambulatorHandler unless not in compute mode");}
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(src_eigvec_index<0)||(src_eigvec_index>=nEigs)){
    errorLaph("Invalid parameters passed to getData in PerambulatorHandler");}
 FileKey fkey(src_time,src_spin);
 RecordKey rkey(snk_spin,snk_time,src_eigvec_index);
 return DHgetPtrSparseGridSummedInsertion->getData(fkey,rkey);
}

Array<std::complex<double>> PerambulatorHandler::getFullData(
                 int snk_time, int snk_spin, int src_time, int src_spin, int nEigsUse) const
{
 if (!isInfoSet()){
    errorLaph("Cannot getData in PerambulatorHandler unless info is set");}
 if (mode==Compute){
    errorLaph("Cannot getData in PerambulatorHandler unless not in compute mode");}
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(nEigsUse<1)||(nEigsUse>nEigs)){
    errorLaph("Invalid parameters passed to getFullData in PerambulatorHandler");}
 FileKey fkey(src_time,src_spin);
 Array<std::complex<double>> result(nEigsUse,nEigsUse);
 for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
    RecordKey rkey(snk_spin,snk_time,srcev_ind);
    const PerambulatorHandler::DataType& buff=DHgetPtr->getData(fkey,rkey);
    for (int n=0;n<nEigsUse;n++){
       result(n,srcev_ind)=buff[n];}
    DHgetPtr->removeData(fkey,rkey);}
 return result;
}

Array<std::complex<double>> PerambulatorHandler::getFullDataSparseGrid(
                 int snk_time, int snk_spin, int src_time, int src_spin, int nEigsUse) const
{
 if (!sgHandler){
    errorLaph("Cannot getFullDataSparseGrid in PerambulatorHandler: no sparse grid specified");}
 if (!isInfoSet()){
    errorLaph("Cannot getData in PerambulatorHandler unless info is set");}
 if (mode==Compute){
    errorLaph("Cannot getData in PerambulatorHandler unless not in compute mode");}
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(nEigsUse<1)||(nEigsUse>nEigs)){
    errorLaph("Invalid parameters passed to getFullData in PerambulatorHandler");}
 FileKey fkey(src_time,src_spin);
 int gridSize = sgHandler->getGrid().getNGridPoints();
 Array<std::complex<double>> result(gridSize,nEigsUse);
 for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
    RecordKey rkey(snk_spin,snk_time,srcev_ind);
    const PerambulatorHandler::DataType& buff=DHgetPtrSparseGrid->getData(fkey,rkey);
    for (int n=0;n<gridSize;n++){
       result(n,srcev_ind)=buff[n];}
    DHgetPtrSparseGrid->removeData(fkey,rkey);}
 return result;
}

Array<std::complex<double>> PerambulatorHandler::getFullDataInsertion(
                 int snk_time, int snk_spin, int src_time, int src_spin, int nEigsUse) const
{
 if (!insHandler){
    errorLaph("Cannot getFullDataInsertion in PerambulatorHandler: no insertion specified");}
 if (!isInfoSet()){
    errorLaph("Cannot getData in PerambulatorHandler unless info is set");}
 if (mode==Compute){
    errorLaph("Cannot getData in PerambulatorHandler unless not in compute mode");}
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(nEigsUse<1)||(nEigsUse>nEigs)){
    errorLaph("Invalid parameters passed to getFullData in PerambulatorHandler");}
 FileKey fkey(src_time,src_spin);
 Array<std::complex<double>> result(nEigsUse,nEigsUse);
 for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
    RecordKey rkey(snk_spin,snk_time,srcev_ind);
    const PerambulatorHandler::DataType& buff=DHgetPtrSummedInsertion->getData(fkey,rkey);
    for (int n=0;n<nEigsUse;n++){
       result(n,srcev_ind)=buff[n];}
    DHgetPtrSummedInsertion->removeData(fkey,rkey);}
 return result;
}

Array<std::complex<double>> PerambulatorHandler::getFullDataSparseGridSummedInsertion(
                 int snk_time, int snk_spin, int src_time, int src_spin, int nEigsUse) const
{
 if (!sginsHandler){
    errorLaph("Cannot getFullDataSparseGridSummedInsertion in PerambulatorHandler: no sparse grid projected insertion specified");}
 if (!isInfoSet()){
    errorLaph("Cannot getData in PerambulatorHandler unless info is set");}
 if (mode==Compute){
    errorLaph("Cannot getData in PerambulatorHandler unless not in compute mode");}
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(nEigsUse<1)||(nEigsUse>nEigs)){
    errorLaph("Invalid parameters passed to getFullData in PerambulatorHandler");}
 FileKey fkey(src_time,src_spin);
 int gridSize = sginsHandler->getGrid().getNGridPoints();
 Array<std::complex<double>> result(gridSize,nEigsUse);
 for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
    RecordKey rkey(snk_spin,snk_time,srcev_ind);
    const PerambulatorHandler::DataType& buff=DHgetPtrSparseGridSummedInsertion->getData(fkey,rkey);
    for (int n=0;n<gridSize;n++){
       result(n,srcev_ind)=buff[n];}
    DHgetPtrSparseGridSummedInsertion->removeData(fkey,rkey);}
 return result;
}

bool PerambulatorHandler::queryData(int snk_time, int snk_spin, int src_time, int src_spin,
                                    int src_eigvec_index) const
{
 if (!isInfoSet()) return false;
 if (mode==Compute) return false;
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(src_eigvec_index<0)||(src_eigvec_index>=nEigs)){
    return false;}
 FileKey fkey(src_time,src_spin);
 RecordKey rkey(snk_spin,snk_time,src_eigvec_index);
 return DHgetPtr->queryData(fkey,rkey);
}

bool PerambulatorHandler::queryDataSparseGrid(int snk_time, int snk_spin, int src_time, int src_spin,
                                    int src_eigvec_index) const
{
 if (!sgHandler){
    errorLaph("Cannot queryDataSparseGrid in PerambulatorHandler: no sparse grid specified");}
 if (!isInfoSet()) return false;
 if (mode==Compute) return false;
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(src_eigvec_index<0)||(src_eigvec_index>=nEigs)){
    return false;}
 FileKey fkey(src_time,src_spin);
 RecordKey rkey(snk_spin,snk_time,src_eigvec_index);
 return DHgetPtrSparseGrid->queryData(fkey,rkey);
}

bool PerambulatorHandler::queryDataInsertion(int snk_time, int snk_spin, int src_time, int src_spin,
                                    int src_eigvec_index) const
{
 if (!insHandler){
    errorLaph("Cannot queryDataInsertion in PerambulatorHandler: no insertion specified");}
 if (!isInfoSet()) return false;
 if (mode==Compute) return false;
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(src_eigvec_index<0)||(src_eigvec_index>=nEigs)){
    return false;}
 FileKey fkey(src_time,src_spin);
 RecordKey rkey(snk_spin,snk_time,src_eigvec_index);
 return DHgetPtrSummedInsertion->queryData(fkey,rkey);
}

bool PerambulatorHandler::queryDataSparseGridSummedInsertion(int snk_time, int snk_spin, int src_time, int src_spin,
                                    int src_eigvec_index) const
{
 if (!sginsHandler){
    errorLaph("Cannot queryDataSparseGridSummedInsertion in PerambulatorHandler: no sparse grid projected insertion specified");}
 if (!isInfoSet()) return false;
 if (mode==Compute) return false;
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(src_eigvec_index<0)||(src_eigvec_index>=nEigs)){
    return false;}
 FileKey fkey(src_time,src_spin);
 RecordKey rkey(snk_spin,snk_time,src_eigvec_index);
 return DHgetPtrSparseGridSummedInsertion->queryData(fkey,rkey);
}


bool PerambulatorHandler::queryFullData(int snk_time, int snk_spin, 
                                        int src_time, int src_spin, int nEigsUse) const
{
 if (!isInfoSet()) return false;
 if (mode==Compute) return false;
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(nEigsUse<1)||(nEigsUse>nEigs)){
    return false;}
 FileKey fkey(src_time,src_spin);
 for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
    RecordKey rkey(snk_spin,snk_time,srcev_ind);
    if (!DHgetPtr->queryData(fkey,rkey)) return false;}
 return true;
}

bool PerambulatorHandler::queryFullDataSparseGrid(int snk_time, int snk_spin,
                                        int src_time, int src_spin, int nEigsUse) const
{
	if (!sgHandler){
		errorLaph("Cannot queryFullDataSparseGrid in PerambulatorHandler: no sparse grid specified");}
	if (!isInfoSet()) return false;
	if (mode==Compute) return false;
	int Textent = uPtr->getTimeExtent();
	int nEigs=getNumberOfLaplacianEigenvectors();
	if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
			||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
			||(nEigsUse<1)||(nEigsUse>nEigs)){
		return false;}
	FileKey fkey(src_time,src_spin);
	int gridSize = sgHandler->getGrid().getNGridPoints();
	for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
		RecordKey rkey(snk_spin,snk_time,srcev_ind);
		if (!DHgetPtrSparseGrid->queryData(fkey,rkey)) return false;}
	return true;
}

bool PerambulatorHandler::queryFullDataInsertion(int snk_time, int snk_spin,
                                        int src_time, int src_spin, int nEigsUse) const
{
 if (!insHandler){
                errorLaph("Cannot queryFullDataInsertion in PerambulatorHandler: no insertion specified");}
 if (!isInfoSet()) return false;
 if (mode==Compute) return false;
 int Textent = uPtr->getTimeExtent();
 int nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
    ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
    ||(nEigsUse<1)||(nEigsUse>nEigs)){
    return false;}
 FileKey fkey(src_time,src_spin);
 for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
    RecordKey rkey(snk_spin,snk_time,srcev_ind);
    if (!DHgetPtrSummedInsertion->queryData(fkey,rkey)) return false;}
 return true;
}

bool PerambulatorHandler::queryFullDataSparseGridSummedInsertion(int snk_time, int snk_spin,
                                        int src_time, int src_spin, int nEigsUse) const
{
        if (!sginsHandler){
                errorLaph("Cannot queryFullDataSparseGrid in PerambulatorHandler: no sparse grid projected insertion specified");}
        if (!isInfoSet()) return false;
        if (mode==Compute) return false;
        int Textent = uPtr->getTimeExtent();
        int nEigs=getNumberOfLaplacianEigenvectors();
        if ((src_time<0)||(src_time>=Textent)||(snk_time<0)||(snk_time>=Textent)
                        ||(src_spin<1)||(src_spin>int(Nspin))||(snk_spin<1)||(snk_spin>int(Nspin))
                        ||(nEigsUse<1)||(nEigsUse>nEigs)){
                return false;}
        FileKey fkey(src_time,src_spin);
        int gridSize = sginsHandler->getGrid().getNGridPoints();
        for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
                RecordKey rkey(snk_spin,snk_time,srcev_ind);
                if (!DHgetPtrSparseGridSummedInsertion->queryData(fkey,rkey)) return false;}
        return true;
}

        // merge data

void PerambulatorHandler::mergeData(const FileListInfo& input_files)
{
 if (!isInfoSet()){
    errorLaph("Cannot mergeData in PerambulatorHandler unless info is set");}
 if (mode!=Merge){
    errorLaph("Cannot mergeData in PerambulatorHandler unless in merge mode");}
 DataGetHandlerMF<PerambulatorHandler,FileKey,RecordKey,DataType> DHget(
           *this,input_files,"Laph--QuarkPeramb","PerambulatorHandlerDataFile");
 set<FileKey> fkeys(DHget.getFileKeys());
 for (set<FileKey>::iterator it=fkeys.begin();it!=fkeys.end();++it){
    DHputPtr->open(*it);
    set<RecordKey> rkeys(DHget.getKeys(*it));
    for (set<RecordKey>::iterator rt=rkeys.begin();rt!=rkeys.end();++rt){
        const DataType& buff=DHget.getData(*it,*rt);
        DHputPtr->putData(*rt,buff);
        DHget.removeData(*it,*rt);}}
}

void PerambulatorHandler::mergeDataSparseGrid(const FileListInfo& input_files)
{
	if (!sgHandler){
		errorLaph("Cannot mergeDataSparseGrid in PerambulatorHandler: no sparse grid specified");}
	if (!isInfoSet()){
		errorLaph("Cannot mergeData in PerambulatorHandler unless info is set");}
	if (mode!=Merge){
		errorLaph("Cannot mergeData in PerambulatorHandler unless in merge mode");}
	DataGetHandlerMF<SparseGridHandler,FileKey,RecordKey,DataType> DHget(
			*sgHandler,input_files,"Laph--SparseGridQuarkPeramb","PerambulatorHandlerSparseGridDataFile");
	set<FileKey> fkeys(DHget.getFileKeys());
	for (set<FileKey>::iterator it=fkeys.begin();it!=fkeys.end();++it){
		DHputPtrSparseGrid->open(*it);
		set<RecordKey> rkeys(DHget.getKeys(*it));
		for (set<RecordKey>::iterator rt=rkeys.begin();rt!=rkeys.end();++rt){
			const DataType& buff=DHget.getData(*it,*rt);
			DHputPtrSparseGrid->putData(*rt,buff);
			DHget.removeData(*it,*rt);}}
}

void PerambulatorHandler::mergeDataSummedInsertion(const FileListInfo& input_files)
{
        if (!insHandler){
                errorLaph("Cannot mergeDataSummedInsertion in PerambulatorHandler: no insertion specified");}
        if (!isInfoSet()){
                errorLaph("Cannot mergeData in PerambulatorHandler unless info is set");}
        if (mode!=Merge){
                errorLaph("Cannot mergeData in PerambulatorHandler unless in merge mode");}
        DataGetHandlerMF<InsertionHandler,FileKey,RecordKey,DataType> DHget(
                        *insHandler,input_files,"Laph--InsertionQuarkPeramb","PerambulatorHandlerDataFile");
        set<FileKey> fkeys(DHget.getFileKeys());
        for (set<FileKey>::iterator it=fkeys.begin();it!=fkeys.end();++it){
                DHputPtrSummedInsertion->open(*it);
                set<RecordKey> rkeys(DHget.getKeys(*it));
                for (set<RecordKey>::iterator rt=rkeys.begin();rt!=rkeys.end();++rt){
                        const DataType& buff=DHget.getData(*it,*rt);
                        DHputPtrSummedInsertion->putData(*rt,buff);
                        DHget.removeData(*it,*rt);}}
}

void PerambulatorHandler::mergeDataSparseGridSummedInsertion(const FileListInfo& input_files)
{
        if (!sginsHandler){
                errorLaph("Cannot mergeDataSparseGridSummedInsertion in PerambulatorHandler: no sparse grid projected insertion specified");}
        if (!isInfoSet()){
                errorLaph("Cannot mergeData in PerambulatorHandler unless info is set");}
        if (mode!=Merge){
                errorLaph("Cannot mergeData in PerambulatorHandler unless in merge mode");}
        DataGetHandlerMF<SparseGridInsertionHandler,FileKey,RecordKey,DataType> DHget(
                        *sginsHandler,input_files,"Laph--SGInsQuarkPeramb","PerambulatorHandlerSparseGridDataFile");
        set<FileKey> fkeys(DHget.getFileKeys());
        for (set<FileKey>::iterator it=fkeys.begin();it!=fkeys.end();++it){
                DHputPtrSparseGridSummedInsertion->open(*it);
                set<RecordKey> rkeys(DHget.getKeys(*it));
                for (set<RecordKey>::iterator rt=rkeys.begin();rt!=rkeys.end();++rt){
                        const DataType& buff=DHget.getData(*it,*rt);
                        DHputPtrSparseGridSummedInsertion->putData(*rt,buff);
                        DHget.removeData(*it,*rt);}}
}
        // check data
   
void PerambulatorHandler::setChecks(const XMLHandler& xmlin)
{
 if (!isInfoSet()){
    errorLaph("Cannot setChecks in PerambulatorHandler unless info is set");}
 if (mode!=Check){
    errorLaph("Cannot setChecks in PerambulatorHandler unless in check mode");}
 XMLHandler xmlrdr(xmlin);
 if (!perambComps.computations.empty()) perambComps.computations.clear();
 uint Textent = uPtr->getTimeExtent();
 uint nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 if (xml_tag_count(xmlrdr,"CheckSet")!=1){
    errorLaph("Cannot setCheckSet in PerambulatorChecker since no <CheckSet> tag");}
 XMLHandler xmlrd(xmlrdr,"CheckSet");
 list<XMLHandler> xmlcs(xmlrd.find("Check"));
 for (list<XMLHandler>::iterator it=xmlcs.begin();it!=xmlcs.end();++it){
    int source_time,tmin,tmax;
    xmlread(*it,"SourceTime",source_time,"LAPH_CHECK_PERAMBULATORS");
    tmin=0; tmax=Textent-1;  
    xmlreadif(*it,"SinkMinimumTime",tmin,"LAPH_PERAMBULATORS");
		xmlreadif(*it,"SinkMaximumTime",tmax,"LAPH_PERAMBULATORS");
		if ((source_time<0)||(source_time>=int(Textent))){
       errorLaph(make_strf("Invalid source time %d",source_time));}
		int nSinkTimes=tmax-tmin+1;
		if ((tmin<=(-int(Textent)))||(tmax>=int(Textent))||(nSinkTimes<=0)||(nSinkTimes>int(Textent))){
			errorLaph(make_strf("Invalid sink tmin and tmax (%d,%d)",tmin,tmax));}
		set<int> srcev_indices;
    if (xml_tag_count(*it,"SourceLaphEigvecIndices")==1){
       vector<int> srcev_inds;
       xmlread(*it,"SourceLaphEigvecIndices",srcev_inds,"LAPH_CHECK_PERAMBULATORS");
       for (int k=0;k<int(srcev_inds.size());++k){
          if ((srcev_inds[k]<0)||(srcev_inds[k]>=int(nEigs))){
             errorLaph(make_strf("Invalid src laph eigvec index %d",srcev_inds[k]));}
       srcev_indices.insert(srcev_inds[k]);}}
    else if (  (xml_tag_count(*it,"SourceLaphEigvecIndexMin")==1)
             &&(xml_tag_count(*it,"SourceLaphEigvecIndexMax")==1)){
       int sevmin=-1, sevmax=-1;
       xmlread(*it,"SourceLaphEigvecIndexMin",sevmin,"LAPH_CHECK_PERAMBULATORS");
       xmlread(*it,"SourceLaphEigvecIndexMax",sevmax,"LAPH_CHECK_PERAMBULATORS");
       if (sevmin<0) sevmin=0;
       if (sevmax>=int(nEigs)) sevmax=nEigs-1;
       if (sevmax<sevmin) sevmax=sevmin;
       for (int ev=sevmin;ev<=sevmax;++ev){
          srcev_indices.insert(ev);}}
    else{
       for (int ev=0;ev<int(nEigs);++ev){
          srcev_indices.insert(ev);}}
    if (srcev_indices.empty()){
       errorLaph("Empty src laph eigvec indices set");}
    perambComps.computations.push_back(PerambComputation(source_time,tmin,tmax,srcev_indices));}
}

   

void PerambulatorHandler::doChecks(const std::string& logfilestub, bool verbose_output)
{
 printLaph("Performing quark perambulator checks");
 int myrank=LayoutInfo::getMyRank();
 int numranks=LayoutInfo::getNumRanks();
 int count=0;
 for (list<PerambComputation>::const_iterator it=perambComps.computations.begin();
      it!=perambComps.computations.end();++it){
    for (uint src_spin=1;src_spin<=Nspin;++src_spin,++count){
       string logfile(logfilestub);
       logfile+=string("_")+make_string(count)+string(".checklog");
       printLaph(make_str("Check results will be sent to ",logfile));
       if (sgHandler) printLaph(make_strf("Check sparse grid results will be sent to %s_SG",logfile));
       if (insHandler) printLaph(make_strf("Check insertion  results will be sent to %s_INS",logfile));
       if (sginsHandler) printLaph(make_strf("Check sparse grid projected insertion results will be sent to ",logfile,"_SGINS"));
       if (count%numranks==myrank){
          doACheck(*it,src_spin,logfile,verbose_output);
		if (sgHandler) doACheckSparseGrid(*it,src_spin,logfile+"_SG",verbose_output);
       		if (insHandler) doACheckInsertion(*it,src_spin,logfile+"_INS",verbose_output);
		if (sginsHandler) doACheckSparseGridSummedInsertion(*it,src_spin,logfile+"_SGINS",verbose_output);}}}
#ifdef ARCH_PARALLEL
 comm_barrier();
#endif
}


   //  DHGetPtr should be in local mode

void PerambulatorHandler::doACheck(const PerambComputation& pcomp, int src_spin,
                                   const std::string& logfile, bool verbose_output)
{
 XMLHandler xmlout("PerambulatorDataCheck");
 xmlout.put_child("SourceTime",make_string(pcomp.src_time));
 xmlout.put_child("SourceSpin",make_string(src_spin));
 string srcevindstr;
 const set<int>& indices(pcomp.src_lapheigvec_indices);
 set<int>::const_iterator vt=indices.begin();
 int rangestart=*vt; int rangestop=*vt+1; ++vt;
 bool keep_going=true;
 while (keep_going){
    if ((vt!=indices.end())&&(rangestop==(*vt))){   // in a started range, increase rangestop
       rangestop++; ++vt;}
    else{                    // end of a range
       if (rangestop==(rangestart+1)){             // only one in range
          srcevindstr+=make_string(rangestart);}
       else{                                       // more than one in range
          srcevindstr+=make_string(rangestart)+"-"+make_string(rangestop-1);}
       if (vt!=indices.end()){
          rangestart=*vt; rangestop=rangestart+1;
          srcevindstr+=" "; ++vt;}
       else{
          keep_going=false;}}}
 xmlout.put_child("SourceLaphEigvecIndices",srcevindstr);

 ofstream fout(logfile);
 fout << "Perambulator check: "<<endl;
 fout << xmlout.output() <<endl;
 uint Textent = uPtr->getTimeExtent();
 uint nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 uint src_time=pcomp.src_time;
 bool flag=true;
 FileKey fkey(src_time,src_spin);
 if (!DHgetPtr->queryFile(fkey)){
    if (verbose_output) fout << " Could not find data for src spin "<<src_spin<<endl;
    flag=false;}
 else{
    for (uint sink_time=0;sink_time<Textent;++sink_time){
    for (uint sinkspin=1;sinkspin<=Nspin;++sinkspin){
    for (set<int>::iterator st=pcomp.src_lapheigvec_indices.begin();
         st!=pcomp.src_lapheigvec_indices.end();++st){
       uint srcev_ind=*st;
       RecordKey rkey(sinkspin,sink_time,srcev_ind);
       if (!DHgetPtr->queryData(fkey,rkey)){
          if (verbose_output) fout << "Could not find record for "<<rkey.output()<<endl;
          flag=false;}
       else{
          const DataType& data=DHgetPtr->getData(fkey,rkey);
          if (data.size()!=nEigs){
             if (verbose_output) fout << "Incorrect data size for record key "<<rkey.output()<<endl;
             flag=false;}
          else{
             bool nanflag=false;
             for (uint k=0;k<data.size();k++){
                if ((std::isnan(data[k].real()))||(std::isnan(data[k].imag()))){
                   nanflag=true; break;}}
             if (nanflag){
                if (verbose_output) fout << "Data contains NaNs for record key "<<rkey.output()<<endl;
                flag=false;}}}}}}}

 if (flag){
    fout << "     Checks OK"<<endl;}
 else{
    fout << "     Some checks FAILED"<<endl;}
 fout.close();
}

void PerambulatorHandler::doACheckSparseGrid(const PerambComputation& pcomp, int src_spin,
                                   const std::string& logfile, bool verbose_output)
{
	if (!sgHandler){
		errorLaph("Cannot doACheckSparseGrid in PerambulatorHandler: no sparse grid specified");}
	XMLHandler xmlout("SparseGridPerambulatorDataCheck");
	xmlout.put_child("SourceTime",make_string(pcomp.src_time));
	xmlout.put_child("SourceSpin",make_string(src_spin));
	string srcevindstr;
	const set<int>& indices(pcomp.src_lapheigvec_indices);
	set<int>::const_iterator vt=indices.begin();
	int rangestart=*vt; int rangestop=*vt+1; ++vt;
	bool keep_going=true;
	while (keep_going){
		if ((vt!=indices.end())&&(rangestop==(*vt))){   // in a started range, increase rangestop
			rangestop++; ++vt;}
		else{                    // end of a range
			if (rangestop==(rangestart+1)){             // only one in range
				srcevindstr+=make_string(rangestart);}
			else{                                       // more than one in range
				srcevindstr+=make_string(rangestart)+"-"+make_string(rangestop-1);}
			if (vt!=indices.end()){
				rangestart=*vt; rangestop=rangestart+1;
				srcevindstr+=" "; ++vt;}
			else{
				keep_going=false;}}}
	xmlout.put_child("SourceLaphEigvecIndices",srcevindstr);

	ofstream fout(logfile);
	fout << "Sparse Grid Perambulator check: "<<endl;
	fout << xmlout.output() <<endl;
	uint Textent = uPtr->getTimeExtent();
	uint nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
	uint src_time=pcomp.src_time;
	int gridSize = sgHandler->getGrid().getNGridPoints();
        int nColor = 3;
	int nTot = nColor*gridSize;	
	bool flag=true;
	FileKey fkey(src_time,src_spin);
	if (!DHgetPtrSparseGrid->queryFile(fkey)){
		if (verbose_output) fout << " Could not find data for src spin "<<src_spin<<endl;
		flag=false;}
	else{
		for (uint sink_time=0;sink_time<Textent;++sink_time){
			for (uint sinkspin=1;sinkspin<=Nspin;++sinkspin){
				for (set<int>::iterator st=pcomp.src_lapheigvec_indices.begin();
						st!=pcomp.src_lapheigvec_indices.end();++st){
					uint srcev_ind=*st;
					RecordKey rkey(sinkspin,sink_time,srcev_ind);
					if (!DHgetPtrSparseGrid->queryData(fkey,rkey)){
						if (verbose_output) fout << "Could not find record for "<<rkey.output()<<endl;
						flag=false;}
					else{
						const DataType& data=DHgetPtrSparseGrid->getData(fkey,rkey);
						if (data.size()!=nTot){
							if (verbose_output) fout << "Incorrect data size for record key "<<rkey.output()<<endl;
							flag=false;}
						else{
							bool nanflag=false;
							for (uint k=0;k<data.size();k++){
								if ((std::isnan(data[k].real()))||(std::isnan(data[k].imag()))){
									nanflag=true; break;}}
							if (nanflag){
								if (verbose_output) fout << "Data contains NaNs for record key "<<rkey.output()<<endl;
								flag=false;}}}}}}}

	if (flag){
		fout << "     Checks OK"<<endl;}
	else{
		fout << "     Some checks FAILED"<<endl;}
	fout.close();
}

void PerambulatorHandler::doACheckInsertion(const PerambComputation& pcomp, int src_spin,
                                   const std::string& logfile, bool verbose_output)
{
 if (!insHandler){
                errorLaph("Cannot doACheckInsertion in PerambulatorHandler: no insertion specified");}
 XMLHandler xmlout("InsertionPerambulatorDataCheck");
 xmlout.put_child("SourceTime",make_string(pcomp.src_time));
 xmlout.put_child("SourceSpin",make_string(src_spin));
 string srcevindstr;
 const set<int>& indices(pcomp.src_lapheigvec_indices);
 set<int>::const_iterator vt=indices.begin();
 int rangestart=*vt; int rangestop=*vt+1; ++vt;
 bool keep_going=true;
 while (keep_going){
    if ((vt!=indices.end())&&(rangestop==(*vt))){   // in a started range, increase rangestop
       rangestop++; ++vt;}
    else{                    // end of a range
       if (rangestop==(rangestart+1)){             // only one in range
          srcevindstr+=make_string(rangestart);}
       else{                                       // more than one in range
          srcevindstr+=make_string(rangestart)+"-"+make_string(rangestop-1);}
       if (vt!=indices.end()){
          rangestart=*vt; rangestop=rangestart+1;
          srcevindstr+=" "; ++vt;}
       else{
          keep_going=false;}}}
 xmlout.put_child("SourceLaphEigvecIndices",srcevindstr);

 ofstream fout(logfile);
 fout << "Insertion Perambulator check: "<<endl;
 fout << xmlout.output() <<endl;
 uint Textent = uPtr->getTimeExtent();
 uint nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 uint src_time=pcomp.src_time;
 bool flag=true;
 FileKey fkey(src_time,src_spin);
 if (!DHgetPtrSummedInsertion->queryFile(fkey)){
    if (verbose_output) fout << " Could not find data for src spin "<<src_spin<<endl;
    flag=false;}
 else{
    for (uint sink_time=0;sink_time<Textent;++sink_time){
    for (uint sinkspin=1;sinkspin<=Nspin;++sinkspin){
    for (set<int>::iterator st=pcomp.src_lapheigvec_indices.begin();
         st!=pcomp.src_lapheigvec_indices.end();++st){
       uint srcev_ind=*st;
       RecordKey rkey(sinkspin,sink_time,srcev_ind);
       if (!DHgetPtrSummedInsertion->queryData(fkey,rkey)){
          if (verbose_output) fout << "Could not find record for "<<rkey.output()<<endl;
          flag=false;}
       else{
          const DataType& data=DHgetPtrSummedInsertion->getData(fkey,rkey);
          if (data.size()!=nEigs){
             if (verbose_output) fout << "Incorrect data size for record key "<<rkey.output()<<endl;
             flag=false;}
          else{
             bool nanflag=false;
             for (uint k=0;k<data.size();k++){
                if ((std::isnan(data[k].real()))||(std::isnan(data[k].imag()))){
                   nanflag=true; break;}}
             if (nanflag){
                if (verbose_output) fout << "Data contains NaNs for record key "<<rkey.output()<<endl;
                flag=false;}}}}}}}
 if (flag){
    fout << "     Checks OK"<<endl;}
 else{
    fout << "     Some checks FAILED"<<endl;}
 fout.close();
}

void PerambulatorHandler::doACheckSparseGridSummedInsertion(const PerambComputation& pcomp, int src_spin,
                                   const std::string& logfile, bool verbose_output)
{
        if (!sginsHandler){
                errorLaph("Cannot doACheckSparseGridSummedInsertion in PerambulatorHandler: no sparse grid projected insertion specified");}
        XMLHandler xmlout("SparseGridSummedInsertionPerambulatorDataCheck");
        xmlout.put_child("SourceTime",make_string(pcomp.src_time));
        xmlout.put_child("SourceSpin",make_string(src_spin));
        string srcevindstr;
        const set<int>& indices(pcomp.src_lapheigvec_indices);
        set<int>::const_iterator vt=indices.begin();
        int rangestart=*vt; int rangestop=*vt+1; ++vt;
        bool keep_going=true;
        while (keep_going){
                if ((vt!=indices.end())&&(rangestop==(*vt))){   // in a started range, increase rangestop
                        rangestop++; ++vt;}
                else{                    // end of a range
                        if (rangestop==(rangestart+1)){             // only one in range
                                srcevindstr+=make_string(rangestart);}
                        else{                                       // more than one in range
                                srcevindstr+=make_string(rangestart)+"-"+make_string(rangestop-1);}
                        if (vt!=indices.end()){
                                rangestart=*vt; rangestop=rangestart+1;
                                srcevindstr+=" "; ++vt;}
                        else{
                                keep_going=false;}}}
        xmlout.put_child("SourceLaphEigvecIndices",srcevindstr);

        ofstream fout(logfile);
        fout << "Sparse Grid Insertion Perambulator check: "<<endl;
        fout << xmlout.output() <<endl;
        uint Textent = uPtr->getTimeExtent();
        uint nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
        uint src_time=pcomp.src_time;
        int gridSize = sginsHandler->getGrid().getNGridPoints();
        int nColor = 3;
        int nTot = nColor*gridSize;
	bool flag=true;
        FileKey fkey(src_time,src_spin);
        if (!DHgetPtrSparseGridSummedInsertion->queryFile(fkey)){
                if (verbose_output) fout << " Could not find data for src spin "<<src_spin<<endl;
                flag=false;}
        else{
                for (uint sink_time=0;sink_time<Textent;++sink_time){
                        for (uint sinkspin=1;sinkspin<=Nspin;++sinkspin){
                                for (set<int>::iterator st=pcomp.src_lapheigvec_indices.begin();
                                                st!=pcomp.src_lapheigvec_indices.end();++st){
                                        uint srcev_ind=*st;
                                        RecordKey rkey(sinkspin,sink_time,srcev_ind);
                                        if (!DHgetPtrSparseGridSummedInsertion->queryData(fkey,rkey)){
                                                if (verbose_output) fout << "Could not find record for "<<rkey.output()<<endl;
                                                flag=false;}
                                        else{
                                                const DataType& data=DHgetPtrSparseGridSummedInsertion->getData(fkey,rkey);
                                                if (data.size()!=nTot){
                                                        if (verbose_output) fout << "Incorrect data size for record key "<<rkey.output()<<endl;
                                                        flag=false;}
                                                else{
                                                        bool nanflag=false;
                                                        for (uint k=0;k<data.size();k++){
                                                                if ((std::isnan(data[k].real()))||(std::isnan(data[k].imag()))){
                                                                        nanflag=true; break;}}
                                                        if (nanflag){
                                                                if (verbose_output) fout << "Data contains NaNs for record key "<<rkey.output()<<endl;
                                                                flag=false;}}}}}}}
if (flag){
                fout << "     Checks OK"<<endl;}
        else{
                fout << "     Some checks FAILED"<<endl;}
        fout.close();
}

// ***************************************************************

  //  static pointers (set to null in default constructor)

unique_ptr<QuarkSmearingHandler> PerambulatorHandler::qSmearHandler;
unique_ptr<GaugeConfigurationHandler> PerambulatorHandler::gaugeHandler;

// ***************************************************************

bool SparseGridHandler::checkHeader(XMLHandler& xmlin, int suffix)
{
 XMLHandler xml_in(xmlin);
 if (xml_tag_count(xml_in,"PerambulatorHandlerSparseGridDataFile")!=1) return false;
 XMLHandler xmlr(xml_in,"PerambulatorHandlerSparseGridDataFile");
 GaugeConfigurationInfo gauge_check(xmlr);
 GluonSmearingInfo gsmear_check(xmlr);
 QuarkSmearingInfo qsmear_check(xmlr);
 uint numspin;
 xmlread(xmlr,"NumSpinComponents", numspin, "SparseGridHandler");
 QuarkActionInfo qaction_check(xmlr);
 RandomSparseGrid grid_check(xmlr);  
 try {
    gInfo.checkEqual(gauge_check);
    gsInfo.checkEqual(gsmear_check);
    qsInfo.checkEqual(qsmear_check); 
    if (numspin!=Nspin){
       throw(std::invalid_argument("Perambulator checkEqual failed...NumSpinComponents mismatch"));}
    qInfo.checkEqual(qaction_check); 
	grid.checkEqual(grid_check); 
 }
 catch(const exception& xp){ return false;}
 return true;
}

void SparseGridHandler::writeHeader(XMLHandler& xmlout, 
                                     const PerambulatorHandler::FileKey& fkey,
                                     int suffix) {
 xmlout.set_root("PerambulatorHandlerSparseGridDataFile");
 XMLHandler xmltmp;
 gInfo.output(xmltmp); xmlout.put_child(xmltmp);
 gsInfo.output(xmltmp); xmlout.put_child(xmltmp);
 qsInfo.output(xmltmp); xmlout.put_child(xmltmp);
 qInfo.output(xmltmp); xmlout.put_child(xmltmp);
 fkey.output(xmltmp); xmlout.put_child(xmltmp);
 xmlout.put_child("NumSpinComponents",make_string(Nspin));
 grid.output(xmltmp); xmlout.put_child(xmltmp); 
}

//------------------------------------------------------------------------------

bool InsertionHandler::checkHeader(XMLHandler& xmlin, int suffix)
{
 XMLHandler xml_in(xmlin);
 if (xml_tag_count(xml_in,"PerambulatorHandlerDataFile")!=1) return false;
 XMLHandler xmlr(xml_in,"PerambulatorHandlerDataFile");
 GaugeConfigurationInfo gauge_check(xmlr);
 GluonSmearingInfo gsmear_check(xmlr);
 QuarkSmearingInfo qsmear_check(xmlr);
 uint numspin;
 xmlread(xmlr,"NumSpinComponents", numspin, "InsertionHandler");
 QuarkActionInfo qaction_check(xmlr);
 InsertionInfo insertion_check(xmlr);  
 try {
    gInfo.checkEqual(gauge_check);
    gsInfo.checkEqual(gsmear_check);
    qsInfo.checkEqual(qsmear_check); 
    if (numspin!=Nspin){
       throw(std::invalid_argument("Perambulator checkEqual failed...NumSpinComponents mismatch"));}
    qInfo.checkEqual(qaction_check); 
	 insInfo.checkEqual(insertion_check); 
 }
 catch(const exception& xp){ return false;}
 return true;
}

void InsertionHandler::writeHeader(XMLHandler& xmlout, 
                                     const PerambulatorHandler::FileKey& fkey,
                                     int suffix) {
 xmlout.set_root("PerambulatorHandlerDataFile");
 XMLHandler xmltmp;
 gInfo.output(xmltmp); xmlout.put_child(xmltmp);
 gsInfo.output(xmltmp); xmlout.put_child(xmltmp);
 qsInfo.output(xmltmp); xmlout.put_child(xmltmp);
 qInfo.output(xmltmp); xmlout.put_child(xmltmp);
 fkey.output(xmltmp); xmlout.put_child(xmltmp);
 xmlout.put_child("NumSpinComponents",make_string(Nspin));
 insInfo.output(xmltmp); xmlout.put_child(xmltmp); 
}

//------------------------------------------------------------------------------

bool SparseGridInsertionHandler::checkHeader(XMLHandler& xmlin, int suffix)
{
 XMLHandler xml_in(xmlin);
 if (xml_tag_count(xml_in,"PerambulatorHandlerSparseGridDataFile")!=1) return false;
 XMLHandler xmlr(xml_in,"PerambulatorHandlerSparseGridDataFile");
 GaugeConfigurationInfo gauge_check(xmlr);
 GluonSmearingInfo gsmear_check(xmlr);
 QuarkSmearingInfo qsmear_check(xmlr);
 uint numspin;
 xmlread(xmlr,"NumSpinComponents", numspin, "SparseGridInsertionHandler");
 QuarkActionInfo qaction_check(xmlr);
 RandomSparseGrid grid_check(xmlr); 
 InsertionInfo insertion_check(xmlr);  
 try {
    gInfo.checkEqual(gauge_check);
    gsInfo.checkEqual(gsmear_check);
    qsInfo.checkEqual(qsmear_check); 
    if (numspin!=Nspin){
       throw(std::invalid_argument("Perambulator checkEqual failed...NumSpinComponents mismatch"));}
    qInfo.checkEqual(qaction_check);
    grid.checkEqual(grid_check);  
	 insInfo.checkEqual(insertion_check); 
 }
 catch(const exception& xp){ return false;}
 return true;
}

void SparseGridInsertionHandler::writeHeader(XMLHandler& xmlout, 
                                     const PerambulatorHandler::FileKey& fkey,
                                     int suffix) {
 xmlout.set_root("PerambulatorHandlerSparseGridDataFile");
 XMLHandler xmltmp;
 gInfo.output(xmltmp); xmlout.put_child(xmltmp);
 gsInfo.output(xmltmp); xmlout.put_child(xmltmp);
 qsInfo.output(xmltmp); xmlout.put_child(xmltmp);
 qInfo.output(xmltmp); xmlout.put_child(xmltmp);
 fkey.output(xmltmp); xmlout.put_child(xmltmp);
 xmlout.put_child("NumSpinComponents",make_string(Nspin));
 grid.output(xmltmp); xmlout.put_child(xmltmp);
 insInfo.output(xmltmp); xmlout.put_child(xmltmp); 
}

//--------------------------------------------------------------------------------
 

void PerambulatorHandler::applyInsertion (LattField& outfield, const LattField& infield)
{
 if (insHandler!=nullptr){
        const std::vector<LattField>& gauge_field=gaugeHandler->getData();
	std::string insert_info = (insHandler->getInsertion()).getInsertionName();
	if(insert_info=="WilsonCloverTwist"){
	    LaphEnv::applyWilsonCloverTwist(outfield, infield, gauge_field, *uPtr, *qactionPtr, insHandler->getInsertion());
	}else if (insert_info=="IsospinMassSplitting"){
            LaphEnv::applyIsospinMassSplitting(outfield, infield);
        }else {
	    printLaph(make_strf("The Name provided in InsertionInfo: '%s' was not yet implemented.",insert_info));
            errorLaph("The Name provided in InsertionInfo was not yet implemented in perambulator_handler.");}
  }else if (sginsHandler!=nullptr){
	const std::vector<LattField>& gauge_field=gaugeHandler->getData();
        std::string insert_info = (sginsHandler->getInsertion()).getInsertionName();
	if(insert_info=="WilsonCloverTwist"){
            LaphEnv::applyWilsonCloverTwist(outfield, infield, gauge_field, *uPtr, *qactionPtr, insHandler->getInsertion());
        }else if (insert_info=="IsospinMassSplitting"){
            LaphEnv::applyIsospinMassSplitting(outfield, infield);
        }else{
            errorLaph("The Name provided in InsertionInfo was not yet implemented in perambulator_handler.");}
   }else {errorLaph("applyInsertion was called although no InsertionHandler is defined ");}
}

}
