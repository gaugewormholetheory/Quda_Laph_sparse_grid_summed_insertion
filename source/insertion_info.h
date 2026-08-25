#ifndef INSERTION_INFO_H
#define INSERTION_INFO_H

#include "xml_handler.h"
#include "quark_action_info.h"

namespace LaphEnv {

// ********************************************************************
// *                                                                  *
// *  Class "InsertionInfo" holds information about                   *
// *  insertions, like the mass_splitting- or TBC-insertion. It also  *
// *  holds the additionally necessary FileListInfo.                  *
// *  The constructor requires an XMLHandler. The XML input must      * 
// *  have the following form:                                        *
// *                                                                  *
// *  <InsertionInfo>                                                 *
// *    <Name>Twist</Name>                                            *
// *    <twist_vector>1 1 1 0</twist_vector>.                         *
// *    <FileListInfo>                                                *
// *      <FileNameStub>ud_peramb_ins_cnfg___CFGNUM__</FileNameStub>  *
// *      <MinFileNumber>0</MinFileNumber>                            *
// *      <MaxFileNumber>3000</MaxFileNumber>                         *
// *    </FileListInfo>                                               *
// *  </InsertionInfo>                                                *
// *                                                                  *
// *  Information about the individual tags now follows:              *
// *                                                                  *
// *    <Name>                                                        *
// *       - simple descriptive name of the insertion type            *
// *       - currently available: WilsonCloverTwist,                  *
// *                              IsospingMassSplitting               *
// *                                                                  *
// *    <twist_vector>                                                *
// *       - necessary for twist-type insertions                      *
// *       - derivative of the twist-angle wrt twist-parameter        *
// *                                                                  *
// ********************************************************************


class InsertionInfo
{

  std::string insertion_name;
  std::vector<int> twist_vector;
  std::set<std::string> types = {"WilsonCloverTwist" , "IsospinMassSplitting"};


 public:
  
    InsertionInfo() = default;
 
    InsertionInfo(XMLHandler& xmlin);

    InsertionInfo(const InsertionInfo& rhs);

    InsertionInfo& operator=(const InsertionInfo& rhs);

    ~InsertionInfo(){}

    void checkEqual(const InsertionInfo& rhs) const;

    bool operator==(const InsertionInfo& rhs) const;

    std::string output(int indent = 0) const;

    void output(XMLHandler& xmlout) const;

    std::vector<int> getTwistVector() const {return twist_vector;}
    
    std::string getInsertionName() const {return insertion_name;}
	
    void setQudaInvertParam(QudaInvertParam& invParam,
                                      const QuarkActionInfo& qactioninfo) const;


 private:

    void setTwistVector(XMLHandler& xmlr, std::vector<int>& twist_vector);
    void set_info(XMLHandler& xmlr);
    bool check_insertion_name(const std::string& ins_name);
    bool check_for_twist(const std::string& ins_name);

};


// **********************************************************
}
#endif
