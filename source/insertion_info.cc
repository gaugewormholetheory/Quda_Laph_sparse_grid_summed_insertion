#include "insertion_info.h"
#include "quda_info.h"
#include <sstream>
#include "laph_stdio.h"

using namespace std;
extern QudaVerbosity getVerbosity();

namespace LaphEnv {
  

// *************************************************************

InsertionInfo::InsertionInfo(XMLHandler& xmlin)
{ 
 XMLHandler xml_in(xmlin);
 set_info(xml_in);
}

void InsertionInfo::set_info(XMLHandler& xml_in)
{
 xml_tag_assert(xml_in,"InsertionInfo");
 XMLHandler xmlr(xml_in,"InsertionInfo");
 xmlread(xmlr,"Name",insertion_name,"InsertionInfo");
 if (!check_insertion_name(insertion_name)) {
  printLaph(make_strf("Unknown Insertion Name in InsertionInfo: '%s' ",insertion_name));
  throw(std::invalid_argument("Unknown Insertion Name"));}
 if(check_for_twist(insertion_name)){
   setTwistVector(xmlr, twist_vector);}
}

void InsertionInfo::setTwistVector(XMLHandler& xmlr, vector<int>& twist_vec)
{
 string vec_str;
 vector<int> temp_vec;
 int val;
 xml_tag_assert(xmlr,"twist_vector");
 xmlread(xmlr,"twist_vector",vec_str,"InsertionInfo");
 istringstream iss(vec_str);
 while (iss >> val) {
 temp_vec.push_back(val);}
 if (temp_vec.size()!=4){
  printLaph("Provided twist vector is not 4d");
  throw(std::invalid_argument("Provided twist vector is not 4d"));
 }
 twist_vec = temp_vec;
}

bool InsertionInfo::check_insertion_name(const std::string& ins_name){
 for (string type : types) {
   if (type==ins_name) return true;
 }
 return false;
}

bool InsertionInfo::check_for_twist(const std::string& ins_name){
 for (int i = 0; i < 5;  ++i){
   if ("Twist"[4-i]!= ins_name.at(ins_name.length()-1-i)) return false;
  }
  return true;
}

  // ************************************************************

     // copy constructor

InsertionInfo::InsertionInfo(const InsertionInfo& rhs)
                : insertion_name(rhs.insertion_name), twist_vector(rhs.twist_vector) {}

						
InsertionInfo& InsertionInfo::operator=(const InsertionInfo& rhs)
{
 insertion_name = rhs.insertion_name;
 if(check_for_twist(rhs.insertion_name)){
   twist_vector = rhs.twist_vector;}
 return *this;
}


void InsertionInfo::checkEqual(const InsertionInfo& rhs) const
{
 if (!((*this)==rhs)){
    printLaph("InsertionInfo checkEqual failed");
    printLaph(make_strf("LHS:\n%s\nRHS:\n%s\n",output(),rhs.output()));
    throw(std::invalid_argument("InsertionInfo checkEqual failed"));}
}



bool InsertionInfo::operator==(const InsertionInfo& rhs) const
{
 return ((insertion_name==rhs.insertion_name)&&(twist_vector==rhs.twist_vector));
}


string InsertionInfo::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}


void InsertionInfo::output(XMLHandler& xmlout) const
{
 xmlout.set_root("InsertionInfo");
 xmlout.put_child("Name",insertion_name);
 xmlout.put_child("twist_vector",make_string(twist_vector));
}

void InsertionInfo::setQudaInvertParam(QudaInvertParam& invParam,
                                      const QuarkActionInfo& qactioninfo) const
{
 invParam=newQudaInvertParam();
 invParam.gamma_basis = QUDA_DIRAC_PAULI_GAMMA_BASIS;
 invParam.dirac_order = QUDA_DIRAC_ORDER;
 invParam.dslash_type = QUDA_WILSON_TBCI_DSLASH;
 invParam.mass_normalization = QUDA_MASS_NORMALIZATION;
 invParam.kappa = qactioninfo.getRValues()[2];
 invParam.mass = qactioninfo.getRValues()[1];
 invParam.input_location = QUDA_CPU_FIELD_LOCATION;
 invParam.output_location = QUDA_CPU_FIELD_LOCATION;

 invParam.cpu_prec = QudaInfo::get_cpu_prec();
 invParam.cuda_prec = QudaInfo::get_cuda_prec();
 
 invParam.solution_type = QUDA_MAT_SOLUTION;
 invParam.solve_type = QUDA_DIRECT_PC_SOLVE;
 invParam.matpc_type = QUDA_MATPC_EVEN_EVEN;
 invParam.dagger = QUDA_DAG_NO;
 
 for (int i=0;i<4;i++) invParam.tbc_ins[i]=twist_vector[i];

 invParam.verbosity = getVerbosity();
 invParam.struct_size = sizeof(invParam);
}
// ***********************************************************
  }
