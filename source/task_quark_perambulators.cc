#include "tasks.h"
#include "perambulator_handler.h"
#include "stop_watch.h"
#include <map>

using namespace std;

namespace LaphEnv {

	
// *****************************************************************************************
// *                                                                                       *
// *   Task to compute the LapH quark perambulators and write to file(s).                  *
// *   The input XML is expected to have the following  form:                              *
// *                                                                                       *
// *      <Task>                                                                           *
// *         <Name>LAPH_QUARK_PERAMBULATORS</Name>                                         *
// *         <QuarkPerambulatorInfo>                                                       *
// *            <GaugeConfigurationInfo>                                                   *
// *               <EnsembleName>CLS_A653</EnsembleName>                                   *
// *               <FileFormat>CERN</FileFormat>                                           *
// *               <ConfigType>WilsonImproved</ConfigType>                                 *
// *               <FileName>/path/A653r000n1</FileName>                                   *
// *               <ConfigNumber>1</ConfigNumber>                                          *
// *               <FermionTimeBC>antiperiodic</FermionTimeBC>                             *
// *            </GaugeConfigurationInfo>                                                  *
// *            <GluonStoutSmearingInfo>                                                   *
// *               <LinkIterations>20</LinkIterations>                                     *
// *               <LinkStapleWeight>0.1</LinkStapleWeight>                                *
// *            </GluonStoutSmearingInfo>                                                  *
// *            <QuarkLaphSmearingInfo>                                                    *
// *               <LaphSigmaCutoff>0.32</LaphSigmaCutoff>                                 *
// *               <NumberLaphEigvecs> 32 </NumberLaphEigvecs>                             *
// *            </QuarkLaphSmearingInfo>                                                   *
// *            <SmearedQuarkFileStub>/path/sq_field</SmearedQuarkFileStub> (see below)    *
// *            <QuarkActionInfo>                                                          *
// *               <Name> WILSON_CLOVER </Name>                                            *
// *               <Flavor> ud </Flavor>                                                   *
// *               <Kappa>0.1365716</Kappa>                                                *
// *               <clovCoeff>2.0668577</clovCoeff>                                        *
// *               <TimeBC> antiperiodic </TimeBC>                                         *
// *             </QuarkActionInfo>                                                        *
// *            <FileListInfo>/path/quark_perambs_1</FileNameStub>                         *
// *               <MinFileNumber>0</MinFileNumber>                                        *
// *               <MaxFileNumber>200</MaxFileNumber>                                      *
// *            </FileListInfo>                                                            *
// *                                                                                       *
// *                                                                                       *
// *   (A RandomSparseGridInfo tag is required if sparse-grid perambs are desired)         *
// *            <RandomSparseGridInfo>                                                     *
// *               <RandomSeed>243525</RandomSeed>                                         *
// *               <GridSpacing>8</GridSpacing>                                            *
// *               <FileListInfo>/path/sparse_grid_quark_perambs_1</FileNameStub>          *
// *                  <MinFileNumber>0</MinFileNumber>                                     *
// *                  <MaxFileNumber>200</MaxFileNumber>                                   *
// *               </FileListInfo>                                                         *
// *            </RandomSparseGridInfo>                                                    *
// *                                                                                       *
// *                                                                                       *
// *   (A InsertionInfo tag is required if insertion perambs are desired)                  *
// *            <InsertionInfo>                                                            *
// *               <Name>WilsonCloverTwist</Name>                                            	   *
// *               <twist_vector>1 1 1 0</twist_vector>                                    *
// *               <FileListInfo>/path/insertion_quark_perambs_1</FileNameStub>            *
// *                  <MinFileNumber>0</MinFileNumber>                                     *
// *                  <MaxFileNumber>200</MaxFileNumber>                                   *
// *               </FileListInfo>                                                         *
// *            </InsertionInfo> 	                                                   *
// *                                                                                       *
// *                   									   *
// *   (If sparse grid projected insertion perambs are desired, both a                     *
// *    RandomSparseGridInfo and a InsertionInfo tag are required. In    	           *
// *    that case one of the two must hold two FileListInfo tags. If both hold two,        *
// *    then the second FileListInfo from RandomSparseGridInfo is used per default)        *
// * 											   *
// *                                                                                       *
// *                     								   *
// *            <InverterInfo>                                                             *
// *               <Name>BICGSTAB</Name>                                                   *
// *               <Tolerance>1.0e-11</Tolerance>                                          *
// *               <MaxIterations>35000</MaxIterations>                                    *
// *            </InverterInfo>                                                            *
// *            <ExtraSolutionChecks/> (optional)                                          *
// *         </QuarkPerambulatorInfo>                                                      *
// *         <ComputationSet>                                                              *
// *            <UseMultiSrcInverter> yes </UseMultiSrcInverter> (optional: yes default)   *
// *            <NumSinksBeforeProject>  4 </NumSinksBeforeProject>                        *
// *            <NumSinksInProjectBatch> 4 </NumSinksInProjectBatch>                       *
// *            <NumEigsInProjectBatch>  4 </NumEigsInProjectBatch>                        *
// *            <Computation>                                                              *
// *               <SourceTime>8</SourceTime>                                              *
// *               <SourceLaphEigvecIndices>12 17 23 24 25</SourceLaphEigvecIndices>       *
// *               <SinkMinimumTime>-12</SinkMinimumTime>                                  *
// *               <SinkMaximumTime>12</SinkMaximumTime>                                   * 
// *             (optional: restrict the output of sink times. Need 0<tmax-tmin<=Text+1)   *
// *             (loop is over tind in [tmin,tmax], but actual time is t=(tind+Text)%Text  *
// *            </Computation>                                                             *
// *            <Computation>                                                              *
// *               <SourceTime>12</SourceTime>                                             *
// *               <SourceLaphEigvecIndexMin>24</SourceLaphEigvecIndexMin>                 *
// *               <SourceLaphEigvecIndexMax>32</SourceLaphEigvecIndexMax>                 *
// *            </Computation>                                                             *
// *         </ComputationSet>                                                             *
// *         <Verbosity>full</Verbosity>  (optional: override default)                     *
// *         <PrintCoefficients/>  (optional)                                              *
// *         <ReportGflops/>  (optional)                                                   *
// *      </Task>                                                                          *
// *                                                                                       *
// *   If the tag <SmearedQuarkFileStub> is set, then the LapH eigenvectors will be read   *
// *   from file.  Otherwise, this task expects that the LapH eigenvectors already         *
// *   reside in HostGlobal; otherwise, an error results.                                  *
// *                                                                                       *
// *   The default verbosity from the main program is used, unless overridden by           *
// *   the <Verbosity> tag above.                                                          *
// *                                                                                       *
// *   The following tags relate to control of the computations:                           *
// *      <UseMultiSrcInverter>  -- yes or no (yes can speed up, uses more memory)         *
// *      <NumSinksBeforeProject>  -- number of inversions to do before projecting         *
// *                                  onto Laph evs (the number to do in multi-src)        *
// *      <NumSinksInProjectBatch> -- while projecting onto Laph evs, the number of        *
// *                                  quark sinks to project at a time as one batch        *
// *                                  (generally should be same as value for previous tag) *
// *      <NumEigsInProjectBatch> -- number of Laph evs to project at a time as one batch  *
// *                                                                                       *
// *   If the <PrintCoefficients/> is present, then the coefficients are printed           *
// *   to standard output.                                                                 *
// *                                                                                       *
// *   Often, all of the perambulators for all source Laph eigenvector indices cannot be   *
// *   done in a single run, and so for safety, the data for different sets of source ev   *
// *   indices will get saved into separate files.  These then need to be combined into    *
// *   a single file for subsequent reading to compute correlators.  The <Merge> task      *
// *   is needed to accomplish this.                                                       *
// *                                                                                       *
// *****************************************************************************************


void doLaphQuarkPerambulators(XMLHandler& xmltask)
{
 if (xml_tag_count(xmltask,"QuarkPerambulatorInfo")!=1){
    errorLaph("Must have one <QuarkPerambulatorInfo> tag");}
 XMLHandler xmlr(xmltask,"QuarkPerambulatorInfo");
 GaugeConfigurationInfo gaugeinfo(xmlr);
 GluonSmearingInfo gSmear(xmlr);
 QuarkSmearingInfo qSmear(xmlr);
 string smeared_quark_filestub;
 xmlreadif(xmlr,"SmearedQuarkFileStub",smeared_quark_filestub,"LAPH_QUARK_PERAMBULATORS");
 QuarkActionInfo quark(xmlr);

 std::map<PerambulatorHandler::PerambulatorType, FileListInfo> Perambulators;
 if (xml_tag_count(xmlr,"Perambulators")!=1){
    errorLaph("Exactly one <Perambulators> list needs to be specified.");}
 XMLHandler xml_per(xmlr,"Perambulators");
 std::list<XMLHandler> entries = xml_per.find_among_children("Perambulator");
 if (entries.empty()){
    errorLaph("The <Perambulators> list must contain at least one <Perambulator>.");}
 for (XMLHandler& entry : entries){
    std::string type_name;
    xmlread(entry,"Type",type_name,"Perambulator");
    PerambulatorHandler::PerambulatorType ptype=PerambulatorHandler::toPerambulatorType(type_name);
    if (ptype==PerambulatorHandler::PERAMBULATOR_INVALID_ENUM){
       errorLaph("Unknown perambulator type: '"+type_name+"' was provided");}
    FileListInfo files_entry(entry);
    if (!Perambulators.emplace(ptype,files_entry).second){
       errorLaph("The provided perambulator types need to be unique.");}}

 using PH = PerambulatorHandler;
 auto has=[&](PH::PerambulatorType t){ return Perambulators.count(t)==1; };

 const bool sg_type     = has(PH::PERAMBULATOR_LEFT_SPARSE_GRID);
 const bool ins_type    = has(PH::PERAMBULATOR_LEFT_INSERTION) || has(PH::PERAMBULATOR_RIGHT_INSERTION)
                       || has(PH::PERAMBULATOR_SUMMED_INSERTION);
 const bool sgins_type  = has(PH::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION)
                       || has(PH::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION);

 // what handler must be provided
 const bool need_sg  = sg_type  || sgins_type;
 const bool need_ins = ins_type || sgins_type;

 if (need_sg && xml_tag_count(xmlr,"RandomSparseGridInfo")!=1){
    errorLaph("Sparse grid perambulator types require exactly one <RandomSparseGridInfo>.");}
 if (need_ins && xml_tag_count(xmlr,"InsertionInfo")!=1){
    errorLaph("Insertion perambulator types require exactly one <InsertionInfo>.");}

 // right-insertion types use a different source than all others
 const bool rescaled_src = has(PH::PERAMBULATOR_RIGHT_INSERTION)
                     || has(PH::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION);
 const bool standard_src = has(PH::PERAMBULATOR_STANDARD) || has(PH::PERAMBULATOR_LEFT_SPARSE_GRID)
                     || has(PH::PERAMBULATOR_LEFT_INSERTION) || has(PH::PERAMBULATOR_SUMMED_INSERTION)
                     || has(PH::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION);
 if (rescaled_src && standard_src){
    errorLaph("Right-insertion types cannot be combined with types based on a different source.");}

 InverterInfo invinfo(xmlr);
 bool upper_spin_only=false;
 if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>0){
    upper_spin_only=true;}
 bool extra_soln_check=false;
 if (xml_tag_count(xmlr,"ExtraSolutionChecks")>0){
    extra_soln_check=true;}
    // change verbosity from default if requested
 Verbosity task_verbosity(getVerbosity());
 if (xml_read_if(xmltask,task_verbosity)){
    setVerbosity(task_verbosity.getQudaValue());}
 bool print_coeffs=false;
 if (xml_tag_count(xmltask,"PrintCoefficients")>0){
    print_coeffs=true;}
 bool report_gflops=false;
 if (xml_tag_count(xmltask,"ReportGflops")>0){
    report_gflops=true;}
 PerambulatorHandler::Mode mode=PerambulatorHandler::Compute;

 printLaph("\n");
 printLaph(" ***********************************************************");
 printLaph(" *                                                         *");
 printLaph(" *   Laph Task: Compute the quark perambulators            *");
 printLaph(" *              and write to file                          *");
 printLaph(" *                                                         *");
 printLaph(" ***********************************************************\n");
 printLaph(make_strf("\n%s\n",gaugeinfo.output()));
 printLaph(make_strf("\n\nGluon Smearing:\n%s\n",gSmear.output()));
 printLaph(make_strf("\n\nQuark Smearing:\n%s\n",qSmear.output()));
 if (!smeared_quark_filestub.empty()){
    printLaph(make_strf("SmearedQuarkFileStub: %s",smeared_quark_filestub));}
 printLaph(make_str("\nQuarkAction:\n",quark.output()));
 printLaph(make_str("\nInverter Info:\n",invinfo.output()));
 if (upper_spin_only){
    printLaph("Only upper spin components used");}
 else{
    printLaph("All spin components used");}
 if (extra_soln_check){
    printLaph("Extra solution checks will be performed");}

 // objects the handler keeps references to: declare BEFORE Q so they are destroyed AFTER it
 std::unique_ptr<RandomSparseGrid> sGridPtr;
 InsertionInfo insertion;
 if (need_sg){
    XMLHandler xmlr_sg(xmlr,"RandomSparseGridInfo");
    sGridPtr.reset(new RandomSparseGrid(xmlr_sg));
    printLaph(make_strf("\nSparse Grid:\n%s\n",sGridPtr->output()));}
 if (need_ins){
    XMLHandler xmlr_ins(xmlr,"InsertionInfo");
    insertion=InsertionInfo(xmlr_ins);
    printLaph(make_strf("\nInsertion:\n%s\n",insertion.output()));}

 // create handler
 PerambulatorHandler Q;
 if (need_sg && need_ins){
    Q.setInfo(gaugeinfo,gSmear,qSmear,quark,*sGridPtr,insertion,Perambulators,
              smeared_quark_filestub,upper_spin_only,mode);}
 else if (need_sg){
    Q.setInfo(gaugeinfo,gSmear,qSmear,quark,*sGridPtr,Perambulators,
              smeared_quark_filestub,upper_spin_only,mode);}
 else if (need_ins){
    Q.setInfo(gaugeinfo,gSmear,qSmear,quark,insertion,Perambulators,
              smeared_quark_filestub,upper_spin_only,mode);}
 else{
    Q.setInfo(gaugeinfo,gSmear,qSmear,quark,Perambulators,
              smeared_quark_filestub,upper_spin_only,mode);}
 
 // read the set of computations (time sources, src eigvec indices)
 // as well as the batching parameters
 XMLHandler xmlcmp(xmltask,"ComputationSet");
 Q.setComputationSet(xmlcmp);

 // set the inverter info
 Q.setInverter(invinfo);
 printLaph("Inverter initialized in PerambulatorHandler");
 Q.outputSuffixMap();

 // now do the computations!
 StopWatch outer; outer.start();
 Q.computePerambulators(extra_soln_check,print_coeffs,report_gflops);
 outer.stop();
 printLaph(make_strf("LAPH_PERAMBULATORS: total time = %g secs",
           outer.getTimeInSeconds()));
 printLaph("LAPH_PERAMBULATORS: ran successfully\n"); 
} 


// *****************************************************************************************
// *                                                                                       *
// *   Task to check the LapH quark perambulators previously computed and written to       *
// *   file.  The input XML is expected to have the following  form:                       *
// *                                                                                       *
// *      <Task>                                                                           *
// *         <Name>LAPH_CHECK_PERAMBULATORS</Name>                                         *
// *         <QuarkPerambulatorInfo>                                                       *
// *            <GaugeConfigurationInfo>                                                   *
// *               <EnsembleName>CLS_A653</EnsembleName>                                   *
// *               <FileFormat>CERN</FileFormat>                                           *
// *               <ConfigType>WilsonImproved</ConfigType>                                 *
// *               <FileName>/path/A653r000n1</FileName>                                   *
// *               <ConfigNumber>1</ConfigNumber>                                          *
// *               <FermionTimeBC>antiperiodic</FermionTimeBC>                             *
// *            </GaugeConfigurationInfo>                                                  *
// *            <GluonStoutSmearingInfo>                                                   *
// *               <LinkIterations>20</LinkIterations>                                     *
// *               <LinkStapleWeight>0.1</LinkStapleWeight>                                *
// *            </GluonStoutSmearingInfo>                                                  *
// *            <QuarkLaphSmearingInfo>                                                    *
// *               <LaphSigmaCutoff>0.32</LaphSigmaCutoff>                                 *
// *               <NumberLaphEigvecs> 32 </NumberLaphEigvecs>                             *
// *            </QuarkLaphSmearingInfo>                                                   *
// *            <SmearedQuarkFileStub>/path/sq_field</SmearedQuarkFileStub>                *
// *            <QuarkActionInfo>                                                          *
// *               <Name> WILSON_CLOVER </Name>                                            *
// *               <Flavor> ud </Flavor>                                                   *
// *               <Kappa>0.1365716</Kappa>                                                *
// *               <clovCoeff>2.0668577</clovCoeff>                                        *
// *               <TimeBC> antiperiodic </TimeBC>                                         *
// *             </QuarkActionInfo>                                                        *
// *            <FileListInfo>/path/quark_perambs_1</FileNameStub>                         *
// *               <MinFileNumber>0</MinFileNumber>                                        *
// *               <MaxFileNumber>200</MaxFileNumber>                                      *
// *            </FileListInfo>                                                            *
// *                                                                                       *
// *                                                                                       *
// *   (A RandomSparseGridInfo tag is required to also check sparse-grid perambs)          *
// *            <RandomSparseGridInfo>                                                     *
// *               <RandomSeed>243525</RandomSeed>                                         *
// *               <GridSpacing>8</GridSpacing>                                            *
// *               <FileListInfo>/path/sparse_grid_quark_perambs_1</FileNameStub>          *
// *                  <MinFileNumber>0</MinFileNumber>                                     *
// *                  <MaxFileNumber>200</MaxFileNumber>                                   *
// *               </FileListInfo>                                                         *
// *            </RandomSparseGridInfo>                                                    *
// *                      								   *
// *   						                                           *
// *   (An InsertionInfo tag is required to also check insertion perambs)                  *
// *            <InsertionInfo>		                                                   *
// *               <Name>Twist</Name>                                      	 	   *
// *               <twist_vector>1 1 1 0</twist_vector>                                    *
// *               <FileListInfo>/path/insertion_quark_perambs_1</FileNameStub>            *
// *                  <MinFileNumber>0</MinFileNumber>                                     *
// *                  <MaxFileNumber>200</MaxFileNumber>                                   *
// *               </FileListInfo>                                                         *
// *            </InsertionInfo>	                                                   *
// *                                                                                       *
// *                 								    	   *
// *   (To check sparse grid projected insertion perambs, both a                           *
// *    RandomSparseGridInfo and a InsertionInfo tag are required. In                      *
// *    that case one of the two must hold two FileListInfo tags)			   *
// *                                                                                       *
// *                     							           *
// *                                                                                       *
// *         </QuarkPerambulatorInfo>                                                      *
// *         <CheckSet>                                                                    *
// *            <Check>                                                                    *
// *               <SourceTime>8</SourceTime>                                              *
// *               <SourceLaphEigvecIndices>12 17 23 24 25</SourceLaphEigvecIndices>       *
// *            </Check>                                                                   *
// *            <Check>                                                                    *
// *               <SourceTime>12</SourceTime>                                             *
// *               <SourceLaphEigvecIndexMin>24</SourceLaphEigvecIndexMin>                 *
// *               <SourceLaphEigvecIndexMax>32</SourceLaphEigvecIndexMax>                 *
// *            </Check>                                                                   *
// *         </CheckSet>                                                                   *
// *         <LogFileStub>/path/logfilestub</LogFileStub>                                  *
// *         <VerboseOutput/> (optional)                                                   *
// *      </Task>                                                                          *
// *                                                                                       *
// *****************************************************************************************

void doLaphCheckPerambulators(XMLHandler& xmltask)
{
 if (xml_tag_count(xmltask,"QuarkPerambulatorInfo")!=1){
    errorLaph("Must have one <QuarkPerambulatorInfo> tag");}
 XMLHandler xmlr(xmltask,"QuarkPerambulatorInfo");
 GaugeConfigurationInfo gaugeinfo(xmlr);
 GluonSmearingInfo gSmear(xmlr);
 QuarkSmearingInfo qSmear(xmlr);
 QuarkActionInfo quark(xmlr);
 
 std::map<PerambulatorHandler::PerambulatorType, FileListInfo> Perambulators;
 bool sparse_grid_type = false;
 bool insertion_type = false;
 bool sparse_grid_insertion_type = false;
 if (xml_tag_count(xmlr,"Perambulators")!=1){
	 errorLaph("Exactly one <Perambulators> list needs to be specified.");}
 XMLHandler xml_per(xmlr,"Perambulators");
 std::list<XMLHandler> entries = xml_per.find_among_children("Perambulator");
 if (entries.empty()){
     errorLaph("The <Perambulators> list must contain at least one <Perambulator> specifying a type.");}
 for (XMLHandler& entry : entries){
	 std::string type_name;
	 xmlread(entry,"Type",type_name,"Perambulator");
     PerambulatorHandler::PerambulatorType perambulator_type
         = PerambulatorHandler::toPerambulatorType(type_name);
     if (perambulator_type == PerambulatorHandler::PERAMBULATOR_INVALID_ENUM){
      	 errorLaph("Unknown perambulator type: '" + type_name + "' was provided");}
   	 FileListInfo files_entry(entry);
   	 if (!Perambulators.emplace(perambulator_type,files_entry).second){
    	 errorLaph("The provided perambulator types need to be unique.");}
 }
 if (Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID)){
	 if (xml_tag_count(xmlr,"RandomSparseGridInfo")!=1)
	 	 errorLaph("If perambulators of sparse grid type were chosen, input must include RandomSparseGridInfo.");
	 sparse_grid_type = true;
 }
 if (Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION)||
 	Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)||
	Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)){
	 if (xml_tag_count(xmlr,"InsertionInfo")!=1)
	 	 errorLaph("If perambulators of insertion type were chosen, input must include InsertionInfo.");
	 insertion_type = true;
 }
 if (Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION)||
	Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION)){
	 if (xml_tag_count(xmlr,"InsertionInfo")!=1)
	 	 errorLaph("If perambulators of sparse grid insertion type were chosen, input must include RandomSparseGridInfo and InsertionInfo.");
	 sparse_grid_insertion_type = true;
 }
 if  ((Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)||
 	Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION))&&
	(Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD)||
	Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID)||
	Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION)||
	Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)||
	Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION))){
	 errorLaph("Only those types of perambulators may be selected whose computations are based on the same source."); 
 }

 InverterInfo invinfo(xmlr);
 bool upper_spin_only=false;
 if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>0){
    upper_spin_only=true;}
 bool verbose_output=false;
 if (xml_tag_count(xmltask,"VerboseOutput")>0){
    verbose_output=true;}
 string logfilestub("perambs_check");
 xmlreadif(xmltask,"LogFileStub",logfilestub,"LAPH_CHECK_PERAMBULATORS");
 PerambulatorHandler::Mode mode=PerambulatorHandler::Check;

 printLaph("\n");
 printLaph(" ***********************************************************");
 printLaph(" *                                                         *");
 printLaph(" *   Laph Task: Check the quark perambulators and          *");
 printLaph(" *              output results to log files                *");
 printLaph(" *                                                         *");
 printLaph(" ***********************************************************\n");
 printLaph(make_strf("\n%s\n",gaugeinfo.output()));
 printLaph(make_strf("\n\nGluon Smearing:\n%s\n",gSmear.output()));
 printLaph(make_strf("\n\nQuark Smearing:\n%s\n",qSmear.output()));
 printLaph(make_str("\nQuarkAction:\n",quark.output()));
 if (upper_spin_only){
    printLaph("Only upper spin components used");}
 else{
    printLaph("All spin components used");}
 printLaph(make_strf("LogFileStub: %s",logfilestub));
 //create handler
 PerambulatorHandler Q;
 InsertionInfo insertion;
 if ((sparse_grid_type)&&(!insertion_type)){
         XMLHandler xmlr_sg(xmlr,"RandomSparseGridInfo");
         RandomSparseGrid sGrid(xmlr_sg);
         printLaph(make_strf("\n\nSparse Grid:\n%s\n",sGrid.output()));
         Q.setInfo(gaugeinfo,gSmear,qSmear,quark,sGrid,Perambulators,
                       "",upper_spin_only,mode);
 }  else if ((!sparse_grid_type)&&(insertion_type)) {
         XMLHandler xmlr_ins(xmlr,"InsertionInfo");
         insertion = InsertionInfo(xmlr);
         printLaph(make_strf("\n\nInsertion:\n%s\n",insertion.output()));
         Q.setInfo(gaugeinfo,gSmear,qSmear,quark,insertion,Perambulators,
                       "",upper_spin_only,mode);
 }  else if (sparse_grid_insertion_type) {
         XMLHandler xmlr_sg(xmlr,"RandomSparseGridInfo");
         RandomSparseGrid sGrid(xmlr_sg);
         printLaph(make_strf("\n\nSparse Grid:\n%s\n",sGrid.output()));
         XMLHandler xmlr_ins(xmlr,"InsertionInfo");
         insertion = InsertionInfo(xmlr_ins);
         printLaph(make_strf("\n\nInsertion:\n%s\n",insertion.output()));
         Q.setInfo(gaugeinfo,gSmear,qSmear,quark,sGrid,insertion,Perambulators,
                       "",upper_spin_only,mode);
 }  else {
         Q.setInfo(gaugeinfo,gSmear,qSmear,quark,Perambulators,"",
                         upper_spin_only,mode);
 }
// read the set of computations (time sources, src eigvec indices)
 // as well as the batching parameters
 XMLHandler xmlchk(xmltask,"CheckSet");
 Q.setChecks(xmlchk);

 // now do the checks!
 StopWatch outer; outer.start();
 Q.doChecks(logfilestub,verbose_output);
 outer.stop();
 printLaph(make_strf("LAPH_CHECK_PERAMBULATORS: total time = %g secs",
           outer.getTimeInSeconds()));
 printLaph("LAPH_CHECK_PERAMBULATORS: ran successfully\n"); 
} 


// *****************************************************************************************
// *                                                                                       *
// *   Task to merge the LapH quark perambulators previously computed and written to       *
// *   file.  The input XML is expected to have the following  form:                       *
// *                                                                                       *
// *      <Task>                                                                           *
// *         <Name>LAPH_MERGE_PERAMBULATORS</Name>                                         *
// *         <QuarkPerambulatorInfo>                                                       *
// *            <GaugeConfigurationInfo>                                                   *
// *               <EnsembleName>CLS_A653</EnsembleName>                                   *
// *               <FileFormat>CERN</FileFormat>                                           *
// *               <ConfigType>WilsonImproved</ConfigType>                                 *
// *               <FileName>/path/A653r000n1</FileName>                                   *
// *               <ConfigNumber>1</ConfigNumber>                                          *
// *               <FermionTimeBC>antiperiodic</FermionTimeBC>                             *
// *            </GaugeConfigurationInfo>                                                  *
// *            <GluonStoutSmearingInfo>                                                   *
// *               <LinkIterations>20</LinkIterations>                                     *
// *               <LinkStapleWeight>0.1</LinkStapleWeight>                                *
// *            </GluonStoutSmearingInfo>                                                  *
// *            <QuarkLaphSmearingInfo>                                                    *
// *               <LaphSigmaCutoff>0.32</LaphSigmaCutoff>                                 *
// *               <NumberLaphEigvecs> 32 </NumberLaphEigvecs>                             *
// *            </QuarkLaphSmearingInfo>                                                   *
// *            <SmearedQuarkFileStub>/path/sq_field</SmearedQuarkFileStub>                *
// *            <QuarkActionInfo>                                                          *
// *               <Name> WILSON_CLOVER </Name>                                            *
// *               <Flavor> ud </Flavor>                                                   *
// *               <Kappa>0.1365716</Kappa>                                                *
// *               <clovCoeff>2.0668577</clovCoeff>                                        *
// *               <TimeBC> antiperiodic </TimeBC>                                         *
// *             </QuarkActionInfo>                                                        *
// *            <FileListInfo>/path/quark_perambs_1</FileNameStub>                         *
// *               <MinFileNumber>0</MinFileNumber>                                        *
// *               <MaxFileNumber>200</MaxFileNumber>                                      *
// *            </FileListInfo>                                                            *
// *                                                                                       *
// *                                                                                       *
// *   (A RandomSparseGridInfo tag is required to also merge sparse-grid perambs)          *
// *            <RandomSparseGridInfo>                                                     *
// *               <RandomSeed>243525</RandomSeed>                                         *
// *               <GridSpacing>8</GridSpacing>                                            *
// *               <FileListInfo>/path/sparse_grid_quark_perambs_1</FileNameStub>          *
// *                  <MinFileNumber>0</MinFileNumber>                                     *
// *                  <MaxFileNumber>200</MaxFileNumber>                                   *
// *               </FileListInfo>                                                         *
// *               <InputFileListInfos>                                                    *
// *                    <FileListInfo>/path/sparse_grid_quark_perambs_1</FileNameStub>     *
// *                            <MinFileNumber>0</MinFileNumber>                           *
// *                            <MaxFileNumber>200</MaxFileNumber>                         *
// *                    </FileListInfo>                                                    *
// *                    <FileListInfo>/path/sparse_quark_perambs_2</FileNameStub>          *
// *                            <MinFileNumber>0</MinFileNumber>                           *
// *                            <MaxFileNumber>200</MaxFileNumber>                         *
// *                    </FileListInfo>                                                    *
// *                    ...                                                                *
// *               </InputFileListInfos>                                                   *
// *            </RandomSparseGridInfo>                                                    *
// *                                                                                       *
// *   (A InsertionInfo tag is required to also merge insertion perambs)                   *
// *            <InsertionInfo>                                                            *
// *               <Name>Twist</Name>                                                      *
// *               <twist_vector>1 1 1 0</twist_vector>                                    *
// *               <FileListInfo>/path/insertion_quark_perambs_1</FileNameStub>            *
// *                  <MinFileNumber>0</MinFileNumber>                                     *
// *                  <MaxFileNumber>200</MaxFileNumber>                                   *
// *               </FileListInfo>                                                         *
// *               <InputFileListInfos>                                                    *
// *            	<FileListInfo>/path/insertion_quark_perambs_1</FileNameStub>       *
// *               		<MinFileNumber>0</MinFileNumber>                           *
// *               		<MaxFileNumber>200</MaxFileNumber>                         *
// *            	</FileListInfo>                                                    *
// *            	<FileListInfo>/path/insertion_perambs_2</FileNameStub>             *
// *               		<MinFileNumber>0</MinFileNumber>                           *
// *               		<MaxFileNumber>200</MaxFileNumber>                         *
// *            	</FileListInfo>                                                    *
// *                    ...                                                                *
// *         	   </InputFileListInfos>                                                   *
// *            </InsertionInfo>                                                           *
// *                                                                                       *
// *                                                                                       *
// *   (To merge sparse grid projected insertion perambs, both a                           *
// *    RandomSparseGridInfo and a InsertionInfo tag are required. In                      *
// *    that case one of the two must hold two FileListInfo tags and one of                *
// *    the two must hold two InputFileListInfo tags)                                      *
// *                                                                                       *
// *                                                                                       *
// *         </QuarkPerambulatorInfo>                                                      *
// *         <InputFileListInfos>                                                          *
// *            <FileListInfo>/path/quark_perambs_1</FileNameStub>                         *
// *               <MinFileNumber>0</MinFileNumber>                                        *
// *               <MaxFileNumber>200</MaxFileNumber>                                      *
// *            </FileListInfo>                                                            *
// *            <FileListInfo>/path/quark_perambs_2</FileNameStub>                         *
// *               <MinFileNumber>0</MinFileNumber>                                        *
// *               <MaxFileNumber>200</MaxFileNumber>                                      *
// *            </FileListInfo>                                                            *
// *                    ...                                                                *
// *         </InputFileListInfos>                                                         * 
// *                                                                                       *
// *											   * 
// *      </Task>                                                                          *
// *                                                                                       *
// *****************************************************************************************

void doLaphMergePerambulators(XMLHandler& xmltask)
{
	if (xml_tag_count(xmltask,"QuarkPerambulatorInfo")!=1){
		errorLaph("Must have one <QuarkPerambulatorInfo> tag");}
	XMLHandler xmlr(xmltask,"QuarkPerambulatorInfo");
	GaugeConfigurationInfo gaugeinfo(xmlr);
	GluonSmearingInfo gSmear(xmlr);
	QuarkSmearingInfo qSmear(xmlr);
	QuarkActionInfo quark(xmlr);

	std::map<PerambulatorHandler::PerambulatorType, FileListInfo> Perambulators;
	std::map<PerambulatorHandler::PerambulatorType, std::list<FileListInfo> > InputPerambulators;
	bool sparse_grid_type = false;
 	bool insertion_type = false;
 	bool sparse_grid_insertion_type = false;
 	if (xml_tag_count(xmlr,"Perambulators")!=1){
	 	errorLaph("Exactly one <Perambulators> list needs to be specified.");}
 	XMLHandler xml_per(xmlr,"Perambulators");
 	std::list<XMLHandler> entries = xml_per.find_among_children("Perambulator");
 	if (entries.empty()){
     	errorLaph("The <Perambulators> list must contain at least one <Perambulator> specifying a type.");}
 	for (XMLHandler& entry : entries){
		std::string type_name;
	 	xmlread(entry,"Type",type_name,"Perambulator");
     	PerambulatorHandler::PerambulatorType perambulator_type
         	= PerambulatorHandler::toPerambulatorType(type_name);
     	if (perambulator_type == PerambulatorHandler::PERAMBULATOR_INVALID_ENUM){
      	 	errorLaph("Unknown perambulator type: '" + type_name + "' was provided");}
   	 	FileListInfo files_entry(entry);
   	 	if (!Perambulators.emplace(perambulator_type,files_entry).second){
    	 	errorLaph("The provided perambulator types need to be unique.");}
		XMLHandler entry_in(entry,"InputFileListInfos");
		std::list<XMLHandler> entry_infs = entry_in.find_among_children("FileListInfo");
		list<FileListInfo> infs;
		for (list<XMLHandler>::iterator it=entry_infs.begin();it!=entry_infs.end();++it)
                 infs.push_back(FileListInfo(*it));
		if (!InputPerambulators.emplace(perambulator_type,infs).second){
    	 	errorLaph("The provided perambulator types need to be unique.");}
 	}
 	if (Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID)){
	 	if (xml_tag_count(xmlr,"RandomSparseGridInfo")!=1)
	 	 	errorLaph("If perambulators of sparse grid type were chosen, input must include RandomSparseGridInfo.");
	 	sparse_grid_type = true;
 	}
 	if (Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION)||
 		Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)||
		Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)){
	 	if (xml_tag_count(xmlr,"InsertionInfo")!=1)
	 	 	errorLaph("If perambulators of insertion type were chosen, input must include InsertionInfo.");
	 	insertion_type = true;
 	}
 	if (Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION)||
		Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION)){
	 	if (xml_tag_count(xmlr,"InsertionInfo")!=1)
	 	 	errorLaph("If perambulators of sparse grid insertion type were chosen, input must include RandomSparseGridInfo and InsertionInfo.");
	 	sparse_grid_insertion_type = true;
 	}
 	if  ((Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION)||
 		Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION))&&
		(Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD)||
		Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID)||
		Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION)||
		Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION)||
		Perambulators.count(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION))){
	 	errorLaph("Only those types of perambulators may be selected whose computations are based on the same source."); 
 	}

 	InverterInfo invinfo(xmlr);
	bool upper_spin_only=false;
	if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>0){
		upper_spin_only=true;}
	PerambulatorHandler::Mode mode=PerambulatorHandler::Merge;

	printLaph("\n");
	printLaph(" ***********************************************************");
	printLaph(" *                                                         *");
	printLaph(" *   Laph Task: Merge the quark perambulators              *");
	printLaph(" *                                                         *");
	printLaph(" ***********************************************************\n");
	printLaph(make_strf("\n%s\n",gaugeinfo.output()));
	printLaph(make_strf("\n\nGluon Smearing:\n%s\n",gSmear.output()));
	printLaph(make_strf("\n\nQuark Smearing:\n%s\n",qSmear.output()));
	printLaph(make_str("\nQuarkAction:\n",quark.output()));
	if (upper_spin_only){
		printLaph("Only upper spin components used");}
	else{
		printLaph("All spin components used");}

	 // create handler
	 PerambulatorHandler Q;
	 InsertionInfo insertion;
	 if ((sparse_grid_type)&&(!insertion_type)){
	         XMLHandler xmlr_sg(xmlr,"RandomSparseGridInfo");
	         RandomSparseGrid sGrid(xmlr_sg);
	         printLaph(make_strf("\n\nSparse Grid:\n%s\n",sGrid.output()));
	         Q.setInfo(gaugeinfo,gSmear,qSmear,quark,sGrid,Perambulators,
	                       "",upper_spin_only,mode);
	 }  else if ((!sparse_grid_type)&&(insertion_type)) {
	         XMLHandler xmlr_ins(xmlr,"InsertionInfo");
	         insertion = InsertionInfo(xmlr_ins);
	         printLaph(make_strf("\n\nInsertion:\n%s\n",insertion.output()));
	         Q.setInfo(gaugeinfo,gSmear,qSmear,quark,insertion,Perambulators,
	                       "",upper_spin_only,mode);
	 }  else if (sparse_grid_insertion_type) {
	     XMLHandler xmlr_sg(xmlr,"RandomSparseGridInfo");
	     RandomSparseGrid sGrid(xmlr_sg);
         printLaph(make_strf("\n\nSparse Grid:\n%s\n",sGrid.output()));
		 XMLHandler xmlr_ins(xmlr,"InsertionInfo");
		 insertion = InsertionInfo(xmlr_ins);
	     printLaph(make_strf("\n\nInsertion:\n%s\n",insertion.output()));
	     Q.setInfo(gaugeinfo,gSmear,qSmear,quark,sGrid,insertion,Perambulators,
	         "",upper_spin_only,mode);
	 }  else {
	         Q.setInfo(gaugeinfo,gSmear,qSmear,quark,Perambulators,"",
        	                 upper_spin_only,mode);
	 }


	//Is this necessary here? 
	// read the set of computations (time sources, src eigvec indices)
	// as well as the batching parameters
	//XMLHandler xmlchk(xmltask,"CheckSet");
	//Q.setChecks(xmlchk);

	StopWatch outer; outer.start();
	auto found_std = InputPerambulators.find(PerambulatorHandler::PerambulatorType::PERAMBULATOR_STANDARD);
	if (found_std != InputPerambulators.end()){	
   		for (auto& fl : found_std->second) Q.mergeData(fl);}
	auto found_sg = InputPerambulators.find(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID);
	if (found_sg != InputPerambulators.end()){	
   		for (auto& fl : found_sg->second) Q.mergeDataSparseGrid(fl);}
	auto found_sins = InputPerambulators.find(PerambulatorHandler::PerambulatorType::PERAMBULATOR_SUMMED_INSERTION);
	if (found_sins != InputPerambulators.end()){	
   		for (auto& fl : found_sins->second) Q.mergeDataSummedInsertion(fl);}
	auto found_lins = InputPerambulators.find(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_INSERTION);
	if (found_lins != InputPerambulators.end()){	
   		for (auto& fl : found_lins->second) Q.mergeDataSummedInsertion(fl);}
	auto found_rins = InputPerambulators.find(PerambulatorHandler::PerambulatorType::PERAMBULATOR_RIGHT_INSERTION);
	if (found_rins != InputPerambulators.end()){	
   		for (auto& fl : found_rins->second) Q.mergeDataSummedInsertion(fl);}
	auto found_sgsins = InputPerambulators.find(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_SUMMED_INSERTION);
	if (found_sgsins != InputPerambulators.end()){	
   		for (auto& fl : found_sgsins->second) Q.mergeDataSparseGridSummedInsertion(fl);}
	auto found_sgrins = InputPerambulators.find(PerambulatorHandler::PerambulatorType::PERAMBULATOR_LEFT_SPARSE_GRID_RIGHT_INSERTION);
	if (found_sgrins != InputPerambulators.end()){	
   		for (auto& fl : found_sgrins->second) Q.mergeDataSparseGridSummedInsertion(fl);}

	outer.stop();
	printLaph(make_strf("LAPH_MERGE_PERAMBULATORS: total time = %g secs",
				outer.getTimeInSeconds()));
	printLaph("LAPH_MERGE_PERAMBULATORS: ran successfully\n"); 
} 


// *****************************************************************************
}
