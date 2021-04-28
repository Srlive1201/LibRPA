#include "input.h"
#include <fstream>
#include<iostream>



using namespace std;
const double PI=3.141592654;
const double TWO_PI=6.283185307;
int NBANDS=0;
int NLOCAL=0;
int NSPIN=0;
int natom=0;
int ncell=0;
matrix wg;
double **ekb;
double efermi;
int n_kpoints=0;
int n_irk_points=0;
int kv_nmp[3]={1,1,1};
Vector3<double> *kvec_c;
Matrix3 latvec;
Matrix3 G;
//vector<ComplexMatrix> wfc_k;	
map<size_t,map<size_t,map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Cs;
map<size_t,map<size_t,map<Vector3_Order<double>,std::shared_ptr<ComplexMatrix>>>> Vq;
map<Vector3_Order<double>,double> irk_wight;
map<int,int> atom_nw;
map<int,int> atom_mu;


void READ_AIMS_BAND(const std::string &file_path)
{
    //cout<<"Begin to read aims-band_out"<<endl;
    ifstream infile;
    infile.open(file_path);
    string ik,is,a,b,c,d;
    infile>>n_kpoints;
    infile>>NSPIN;
    infile>>NBANDS;
    infile>>NLOCAL;
    infile>>efermi;
    efermi*=2;
    efermi=0.2;
    //cout<<n_kpoints<<NSPIN<<NLOCAL<<endl;
    wg.create(n_kpoints,NLOCAL);

    ekb=new double*[n_kpoints];
    for(int ik=0;ik!=n_kpoints;ik++)
    {
        ekb[ik]=new double[NLOCAL];
    }

    for(int ks=0;ks!=n_kpoints;ks++)
    {
        infile>>ik>>is;
        //cout<<ik<<is<<endl;
        int k_index=stoi(ik)-1;
        int s_index=stoi(is)-1;
        for(int i=0;i!=NLOCAL;i++)
        {
            infile>>a>>b>>c>>d;
            //cout<<a<<b<<c<<d<<endl;
            wg(k_index,i)=stod(b)/n_kpoints; //different with abacus!
            ekb[k_index][i]=stod(c)*2;
        }

    }
    cout<<"efermi :"<<efermi<<endl;
    cout<<"Success read aims-band_out"<<endl;
} 

void READ_AIMS_STRU(const std::string &file_path)
{
    cout<<"Begin to read aims stru"<<endl;
    ifstream infile;
    string x,y,z;
    infile.open(file_path);
    infile >>x>>y>>z;
    latvec.e11=stod(x);latvec.e12=stod(y);latvec.e13=stod(z);
    infile >>x>>y>>z;
    latvec.e21=stod(x);latvec.e22=stod(y);latvec.e23=stod(z);
    infile >>x>>y>>z;
    latvec.e31=stod(x);latvec.e32=stod(y);latvec.e33=stod(z);
    latvec/=1.889726127;
    latvec.print();

    infile >>x>>y>>z;
    G.e11=stod(x);G.e12=stod(y);G.e13=stod(z);
    infile >>x>>y>>z;
    G.e21=stod(x);G.e22=stod(y);G.e23=stod(z);
    infile >>x>>y>>z;
    G.e31=stod(x);G.e32=stod(y);G.e33=stod(z);
    
    G/=TWO_PI;
    G*=1.889726127;
    G.print();
    Matrix3 latG=latvec*G;
    cout<<" lat * G"<<endl;
    latG.print();
    infile >>x>>y>>z;
    kv_nmp[0]=stoi(x);
    kv_nmp[1]=stoi(y);
    kv_nmp[2]=stoi(z);
    kvec_c=new Vector3<double>[n_kpoints];
    for(int i=0;i!=n_kpoints;i++)
    {
        infile>>x>>y>>z;
        kvec_c[i]={stod(x),stod(y),stod(z)};
        kvec_c[i]*=(1.889726127/TWO_PI);
        cout<<kvec_c[i]<<endl;
    }
}

void READ_AIMS_EIGENVECTOR(const std::string &file_path, vector<ComplexMatrix> &wfc_k)
{
    cout<<"Begin to read aims eigenvecor"<<endl;
    ifstream infile;
    infile.open(file_path);
    string rvalue, ivalue;
    assert(NSPIN==1);
    cout<<" n_kpoints:" <<n_kpoints<<endl;
    wfc_k.resize(n_kpoints);
    cout<<"wfc_k size: "<<wfc_k.size()<<endl;
    for(int ik=0;ik!=n_kpoints;ik++)
        wfc_k[ik].create(NBANDS,NLOCAL);
    cout<<"  create wfc_k"<<endl;
    for(int iw=0;iw!=NLOCAL;iw++)
        for(int ib=0;ib!=NBANDS;ib++)
            for(int is=0;is!=NSPIN;is++)
                for(int ik=0;ik!=n_kpoints;ik++)
                {
                    //cout<<iw<<ib<<is<<ik;
                    infile>>rvalue>>ivalue;
                    //cout<<rvalue<<ivalue<<endl;
                    wfc_k[ik](ib,iw)=complex<double>(stod(rvalue),stod(ivalue));
                }
    //cout<<" EOF?  "<<infile.eof()<<endl;
}

void READ_AIMS_Cs(const std::string &file_path)
{
    //map<size_t,map<size_t,map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Cs_m;
    //cout<<"Begin to read aims Cs"<<endl;
    string natom_s,ncell_s,ia1_s,ia2_s,ic_1,ic_2,ic_3,i_s,j_s,mu_s,Cs_ele;
    ifstream infile;
    infile.open(file_path);
    infile>>natom_s>>ncell_s;
    natom=stoi(natom_s);
    ncell=stoi(ncell_s);
    cout<<"  Natom  Ncell  "<<natom<<ncell<<endl;
    for(int loop=0;loop!=natom*natom*ncell;loop++)
    {
        infile>>ia1_s>>ia2_s>>ic_1>>ic_2>>ic_3>>i_s;
        infile>>j_s>>mu_s;
        //cout<<ic_1<<mu_s<<endl;
        int ia1=stoi(ia1_s)-1;
        int ia2=stoi(ia2_s)-1;
        int ic1=stoi(ic_1);
        int ic2=stoi(ic_2);
        int ic3=stoi(ic_3);
        int n_i=stoi(i_s);
        int n_j=stoi(j_s);
        int n_mu=stoi(mu_s);
        
        atom_nw.insert(pair<int,int>(ia1,n_i));
        atom_mu.insert(pair<int,int>(ia1,n_mu));
        Vector3_Order<int> box(ic1,ic2,ic3);
        //cout<< ia1<<ia2<<box<<endl;

        
        shared_ptr<matrix> cs_ptr=make_shared<matrix>();
        cs_ptr->create(n_i*n_j,n_mu);
        //cout<<cs_ptr->nr<<cs_ptr->nc<<endl;
        
        for(int i=0;i!=n_i;i++)
            for(int j=0;j!=n_j;j++)
                for(int mu=0;mu!=n_mu;mu++)
                {
                    infile>>Cs_ele;
                    (*cs_ptr)(i*n_j+j,mu)=stod(Cs_ele);
                    //(*cs_ptr)(j*n_i+i,mu)=stod(Cs_ele);
                }
        //if(abs(stoi(ic_3))>=10) continue;
        Cs[ia1][ia2][box]=cs_ptr;
    }
    cout<<" atom_mu.size= "<<atom_mu.size()<<endl;
    cout<<"  Cs_dim   "<<Cs.size()<<"    "<<Cs[0].size()<<"   "<<Cs[0][0].size()<<endl;
    cout<<"FINISH to read aims Cs"<<endl;
}

void READ_AIMS_Vq(const std::string &file_path)
{
    cout<<"Begin to read aims vq"<<endl;
    ifstream infile;
    infile.open(file_path);
    string iat_1, iat_2, mu_s, nu_s, q1,q2,q3, vq_r,vq_i;
    for(int i=0;i!=natom*natom;i++)
    {
        infile>> iat_1>>iat_2>>mu_s>>nu_s;
        infile>> q1>>q2>>q3;
        int ia1=stoi(iat_1)-1;
        int ia2=stoi(iat_2)-1;
        int mu=stoi(mu_s);
        int nu=stoi(nu_s);
        double qx=stod(q1);
        double qy=stod(q2);
        double qz=stod(q3);
        Vector3_Order<double> qvec(qx,qy,qz);
        shared_ptr<ComplexMatrix>  vq_ptr=make_shared<ComplexMatrix>();
        vq_ptr->create(mu,nu);
        for(int i_nu=0;i_nu!=nu;i_nu++)
            for(int i_mu=0;i_mu!=mu;i_mu++)
            {
                infile>>vq_r>>vq_i;
                (*vq_ptr)(i_mu,i_nu)=complex<double>(stod(vq_r),stod(vq_i));
            }
        Vq[ia1][ia2][qvec]=vq_ptr;
    }
    cout<<"  Vq_dim   "<<Vq.size()<<"    "<<Vq[0].size()<<"   "<<Vq[0][0].size()<<endl;
    cout<<"Finish read aims vq"<<endl;
}

void READ_AIMS_Vq_real(const std::string &file_path)
{
    cout<<"Begin to read aims vq_real"<<endl;
    ifstream infile;
    infile.open(file_path);
    string c_size,iat_1, iat_2, mu_s, nu_s, q1,q2,q3, vq_r,vq_i,q_num,q_weight;
    infile >> n_irk_points;
    for(int ik=0;ik!=n_irk_points;ik++)
    {
        infile>> c_size>>mu_s>>nu_s;
        infile>> q_num >>q_weight;
        int mu=stoi(mu_s);
        int nu=stoi(nu_s);
        
        irk_wight.insert(pair<Vector3_Order<double>,double>(kvec_c[stoi(q_num)-1],stod(q_weight)));
        shared_ptr<ComplexMatrix>  vq_full_ptr=make_shared<ComplexMatrix>();
        vq_full_ptr->create(mu,nu);
        for(int i_nu=0;i_nu!=nu;i_nu++)
            for(int i_mu=0;i_mu!=mu;i_mu++)
            {
                infile>>vq_r>>vq_i;
                (*vq_full_ptr)(i_mu,i_nu)=complex<double>(stod(vq_r),stod(vq_i));
            }
        //Vq[0][0][Vector3_Order<double>{0.0,0.0,0.0}]=vq_ptr;
        ComplexMatrix vq_full=*vq_full_ptr;
        for(int I=0;I!=atom_mu.size();I++)
            for(int J=0;J!=atom_mu.size();J++)
            {
                shared_ptr<ComplexMatrix>  vq_ptr=make_shared<ComplexMatrix>();
                vq_ptr->create(atom_mu[I],atom_mu[J]);
                for(int i_mu=0;i_mu!=atom_mu[I];i_mu++)
                    for(int i_nu=0;i_nu!=atom_mu[J];i_nu++)
                    {
                        (*vq_ptr)(i_mu,i_nu)=vq_full(atom_mu_loc2glo(I,i_mu),atom_mu_loc2glo(J,i_nu));
                    }
                Vector3_Order<double> qvec(kvec_c[stoi(q_num)-1]);
                cout<<q_num<<qvec<<endl;
                Vq[J][I][qvec]=vq_ptr;
            }
    }
    cout<<"  Vq_dim   "<<Vq.size()<<"    "<<Vq[0].size()<<"   "<<Vq[0][0].size()<<endl;
    for(auto &irk:irk_wight)
        cout<<" irk_vec and weight: "<<irk.first<<"  "<<irk.second<<endl;
    cout<<"Finish read aims vq"<<endl;
}
int atom_iw_loc2glo(const int &atom_index,const int & iw_lcoal)
{
    int nb=0;
    for(int ia=0;ia!=atom_index;ia++)
        nb+=atom_nw[ia]; 
    return iw_lcoal+nb;
}

int atom_mu_loc2glo(const int &atom_index,const int & mu_lcoal)
{
    int nb=0;
    for(int ia=0;ia!=atom_index;ia++)
        nb+=atom_mu[ia]; 
    return mu_lcoal+nb;
}
