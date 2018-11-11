# include <stdio.h>
# include <sys/types.h>
# include <math.h>
# include <stdlib.h>
# include <time.h>
# include <fstream>
# include <iostream>
# define maxNp 100
# define maxN 50
# define SQ(x) ((x)*(x))
# define CUB(x) ((x)*(x)*(x))
using namespace std;
void distance_ads(double,double,double,int,double,double,double,double *,double *,double *,double *);
void distance(int,int,double,double,double,double *,double *,double *,double *);
void overlapvolume(double,double,double,double,double,double,double,double,double ,double *,double *,double *,double *);
void overlapvolume_surface(double,double,double,double ,double *,double *,double *,double *,double *,double *);
void overlap(double,double,double,double *);
void overlap_prime(double,double,double,double *);
void overlap_surface(double,double,double,double,double *);
void overlap_surface_prime(double,double,double,double,double *);
void normal(double *);
typedef struct 
{
	double cum;
	int IPi;			 //Maximum number of free linkers
	int IPj;		     //Maximum number of neighbors
	int type;
}Rates;
Rates R[maxNp];
typedef struct 
{
	double xl,yl,zl;
	int nl;

}Linkers;
Linkers L[maxN];
typedef struct 
{	
	int nn_index; 			              //index of neighbors
    int nb;
	int nb_AB,nb_BA;				    //number of bridges
	int n_index;
	double aoff,aon,aon_AB,aon_BA,aoff_AB,aoff_BA,aoff_L,aon_L;
	double poff,pon,poff_AB,pon_AB,poff_BA,pon_BA,poff_L,pon_L;
    double d;
    //double overlap_i_prime,omega_ij;
    double e_ij;
    double e_ij_prime;
    double omega_ij;
    double omega_ij_prime;
    int nl;
    double nl_d,nb_d;
}Pairs;
Pairs Pn[maxNp];
typedef struct 
{
	double x,y,z;
    int nf;
	int nf_A,nf_B,maxn,maxnn;			//Maximum number of free linkers, Maximum number of neighbors
    
    double nf_A_old,nf_B_old;
    double nf_A_new,nf_B_new;
    double nbS_AC;
    double nbS_d;
    //int nbS_AB,nbS_BA;
    int Sc;                             //index of the neighbor particles of the surface  
    double Sd;                          //Distance between the neighbor and surface
	Linkers L[maxN];
	Pairs Pn[maxNp];
    double omega_i;
    int type;
    double Fx,Fy,Fz;
    double aon_AB_S,aon_BA_S,aoff_AB_S,aoff_BA_S;
    double pon_AB_S,pon_BA_S,poff_AB_S,poff_BA_S;
    double e_is;
    double e_is_prime;
    double omega_is;
    double omega_is_prime;
}Particles;
Particles P[maxNp];

double kon,koff,kon_L,koff_L;


int main()
{   
    int z=12;
    int Np=z+1;
	
    int NA=200;
    int NB=200;
    int N=NA+NB;
    int NsC=400;
    int Ns=NsC;
    //int NsB=100;
    double Rc=20.0;
	double L=1.0;
    double R2s=1.0;
	double Lx=150.0;
	double Ly=150.0;
	double Lz=150.0;
    //double LL=15.0;
    //double LL1=0.0;
    double beta=1.0;
    double rho0=1.0;
	double Delta_tB=0.001;
    double t_tot=1.0;
    double kon0=100000.0;
    double DG0=-10.0;
    double DG0_star=beta*DG0-log(1/rho0*CUB(L));
    double DG0_L=-10.0;
    double DG0_L_star=beta*DG0_L-log(1/rho0*CUB(L));
    
    int nfS_C;
    double nfS_C_old,nfS_C_new;
    
	double dx,dy,dz,d,X,Y,Z;
	double t,t_bar,tao;
    
	double a_on_sum,a_off_sum,a_tot;
    double a_on_AB_sum, a_on_BA_sum, a_off_AB_sum, a_off_BA_sum, a_on_sum_L, a_off_sum_L;
    double a_on_AB_sum_S, a_on_BA_sum_S, a_off_AB_sum_S, a_off_BA_sum_S; 
    int n_deposited_particles, n_checked_particles;
    
	int count,count_nn,count_n;
	//int seed =time(NULL)*getpid();
    int seed = 12;
    srand (seed);
    //cout << "seed=" << seed << endl;
	
	int i,j,k,E,q;
	bool flag;
    //ofstream myfile;
    //ofstream myfile1;
    ofstream myfile2;
    //myfile.open ("Particle_gilles.txt");
    //myfile1.open ("Particle_gilles.dump");
    myfile2.open ("Phase_diagram_b_l.dump");
    
    for (i=0; i<Np; i++){
        P[i].nf_A=NA;
        P[i].nf_B=NB;
    }
    nfS_C=NsC;           //number of free surface linker (initiation) A type 
    //nfS_B=NsB;         //number of free surface linker (initiation) B type
    

/***********************************Particle position *******************************************************************/    
    n_deposited_particles=0;
    for (i=0; i<1; i++)
    {
        P[i].x =(rand()/(double)RAND_MAX)*Lx;  //Lx
        P[i].y =(rand()/(double)RAND_MAX)*Ly;
        P[i].z =(rand()/(double)RAND_MAX)*Lz;
    }
    
    n_deposited_particles=1;

    while (n_deposited_particles < Np)
    {
        X =(rand()/(double)RAND_MAX)*Lx;
        Y =(rand()/(double)RAND_MAX)*Ly;
        Z =(rand()/(double)RAND_MAX)*Lz;
        
        n_checked_particles = 0;
        for (i=0; i<n_deposited_particles; i++)
        {
            distance_ads(X,Y,Z,i,Lx,Ly,Lz,&dx,&dy,&dz,&d);

            if (d > 2*Rc) {
                n_checked_particles = n_checked_particles + 1;
            }
            else if (d <= 2*Rc) break;
        }
        if (n_checked_particles == n_deposited_particles){
            
            P[n_deposited_particles].x=X;
            P[n_deposited_particles].y=Y;
            P[n_deposited_particles].z=Z;
            
            n_deposited_particles = n_deposited_particles + 1;
        } 
    }
//    P[0].x=0.0;     P[0].y=0.0;     P[0].z=0.0;
//    P[1].x=11.5;    P[1].y=0.0;     P[1].z=0.0;
//    P[2].x=0.0;     P[2].y=11.5;    P[2].z=0.0;
//    P[3].x=-11.5;   P[3].y=0.0;     P[3].z=0.0;
//    P[4].x=0.0;     P[4].y=-11.5;   P[4].z=0.0;
//    P[5].x=0.0;     P[5].y=-11.5;   P[5].z=0.0;
//    P[6].x=0.0;     P[6].y=-11.5;   P[6].z=0.0;
    
    
    
    double Sx=Lx/2.0;
    double Sy=Ly/2.0;
    double Sz=0.0;
    //exit(0);
    //cout << "n_deposited_particles=" << n_deposited_particles << endl;
    //myfile << n_deposited_particles <<  endl;
    //myfile << "Atoms. Timestep: 0" <<  endl;
    for(i=0; i<n_deposited_particles; i++)
    {
        //cout << " " << P[i].x << " " << P[i].y << " " << P[i].z <<  endl;
        //printf("%lf %lf %lf\n",P[i].x,P[i].y,P[i].z);
        //myfile << "N" << " " << P[i].x << " " << P[i].y << " " << P[i].z <<  endl;    
    }
/****************************************************neighbor calculations*********************************************************************/
    t=0.0;
	//while (t < t_tot)  //For Brownian dynamics; for the time being it is commented out
	{
        int ns=0;
        for(i=0; i<Np; i++)
        {
            count_n=0;
            for (j=i+1; j<Np; j++)
            {
                distance(i,j,Lx,Ly,Lz,&dx,&dy,&dz,&d);
                d=2*Rc+1.1; //11.1;
                
	 		    if (d < 2*(Rc+L) )                         //Rc=5 && L=1.0
	 		    {
                    count_nn=P[i].maxnn;                     //counter
                    P[i].Pn[count_nn].nn_index=j;		     //index of neighbors
                    P[i].Pn[j].d=d; 
                    P[i].maxnn=P[i].maxnn+1;			     //Maximum number of neighbors
                                                               
                    
                    count_nn=P[j].maxnn;                     //counter
                    P[j].Pn[count_nn].nn_index=i;		     //index of neighbors
                    P[j].Pn[i].d=d; 
                    P[j].maxnn=P[j].maxnn+1;			     //Maximum number of neighbors
                    
                    
/********************************************* To calculate neigbors only one time *************************************************************************/              
                    P[i].Pn[count_n].n_index=j;				 //index of neighbors without repeat
                    count_n=count_n+1;                       //counter
                    P[i].maxn=P[i].maxn+1;			         //Maximum number of neighbors without repeat
/***********************************************************************************************************************************************************/
                }
            }
            
        }
        //cout << P[6].maxnn << endl;
        
        //exit(0);
/**********************************************************************************************************************************************/
        double e=2.71828;
        double vol=(4.0/3.0)*M_PI*CUB(Rc);
        //double rhog=0.01/vol; //2.0;
        double rhog=2.1e-6/e; //1.1e-05; //7.0e-04; //0.00001; // ; 4.36006E-06 4.1e-06
        //cout << log(rhog)*pow((L/2.0),3)+1 << " " << log(0.01)*pow((L/2.0),3)+1 <<endl;
        //exit(0);
        //double rhog=0.01;
        double delta=L; // for delta=L/2,L,2L
        double target=(log(rhog*CUB(delta/2.0))+1); // /(z/2.0);
        //cout << target << " " << rhog << " " << (log(0.01)*CUB(delta/2.0)+1)/(z/2.0) << endl;
        //exit(0);
/***********************************  overlap calculation *******************************************/
        double R1  =Rc+L;
        double R2  =Rc;
        double D1  =R1+R1;
        double D1_ =R1-R1;
        double D2  =R2+R2;
        double D2_ =R2-R2;
        double D12 =R1+R2;
        double D12_=R1-R2;
        double ov, ov1, ov2, ov3, ov_prime, ov1_prime, ov2_prime, ov3_prime;
        double e_ij,w_ij,e_ij_prime,w_ij_prime;
        //double overlap_sum_i,overlap_i,overlap_ij,overlap_i_prime,overlap_ij_prime;
        double e_ij_sum;
        double omega_0=4*M_PI*SQ(Rc)*L;
        double r;
        for (int counter_r=1; counter_r<2; counter_r++)
        {
        r=2*Rc+1.0; //11.0; 
        for (int counter_DG0=0; counter_DG0<19; counter_DG0++)
        {
            //r=7.0-counter_r*0.04;
            DG0_star=-20.0+counter_DG0;
            DG0_L_star=DG0_star; //log(0.15); //+log(0.15); //+log((314.159/5.89782)*(1/rho0));
            
            //DG0_star=-15.75;
            //DG0_L_star=DG0_star; //-7.4; //DG0_star+1;
            for (int counter_N=0; counter_N<2000; counter_N=counter_N+10) //200   +2
            {
                //N=200;
                N=2+counter_N;
                NA=N/2.0;
                NB=N/2.0;
                cout << NA << " " << NB << endl;
                //exit(0);
               /*****colloid-colloid overlap function****/
                for (i=0; i<Np; i++)
                {   
                    e_ij_sum=0.0;
                    for (j=0; j<P[i].maxnn; j++)
                    {
                        int index = P[i].Pn[j].nn_index;
                        //double r=P[i].Pn[index].d;
                        overlapvolume(r,Rc,L,D1,D1_,D2,D2_,D12,D12_,&e_ij,&w_ij,&e_ij_prime,&w_ij_prime);
                        e_ij_sum=e_ij_sum+e_ij;
                
                        P[i].Pn[index].e_ij = e_ij;
                        P[i].Pn[index].e_ij_prime = e_ij_prime;
                        P[i].Pn[index].omega_ij= w_ij;
                        //P[i].Pn[index].omega_ij = omega_0*0.1;
                        P[i].Pn[index].omega_ij_prime = w_ij_prime; 
                    }
                    P[i].omega_i=omega_0-e_ij_sum;
                }
                
                cout << P[0].omega_i << " " << omega_0 << " " << P[0].Pn[2].omega_ij/P[2].omega_i << endl;
                //exit(0);

/***************************** Fixpoint iteration ****************************************/
    
            for (i=0; i<Np; i++)
            {   
                P[i].nf_A_old=rand()% NA+1;
                P[i].nf_B_old=rand()% NB+1;
                P[i].nf_A_new=P[i].nf_A_old;
                P[i].nf_B_new=P[i].nf_B_old;
                //cout << P[i].nf_A_old << endl;
            }
            //nfS_C_old=rand()% NsC+1;
            //nfS_C_new=nfS_C_old;
            
            //nfS_B_new=nfS_B_old;
            //cout << " Random " << P[0].nf_A_old << " " << P[0].nf_B_old << endl;
            //exit(0);
            double sum,sum_B;
            flag=true;
            double Np_roots;
            int index;
            int Niteration=0;
            double error,error_B,error_C;
            double Free_energy_regularization;
            double betaDG_ij;
            double betaDG_ii_L;
            while (flag==true && Niteration < 10000000)
            {
                Np_roots=0.0;
                for (i=0; i<Np; i++)
                {
                    sum=0.0;
                    sum_B=0.0;
                    for (j=0; j<P[i].maxnn; j++)
                    {
                        index=P[i].Pn[j].nn_index;
                       
                        
                        double PF = P[i].Pn[index].omega_ij/(P[i].omega_i*P[index].omega_i);
                        
                        double DG_ij=DG0_star/beta-log(PF)/beta;
                        betaDG_ij=beta*DG_ij;
                    
                        sum=sum+P[index].nf_B_old*exp(-betaDG_ij);
                        sum_B=sum_B+P[index].nf_A_old*exp(-betaDG_ij);
                    }
                    double PF_L = 1/(P[i].omega_i*rho0);
                    double DG_ii_L = DG0_L_star/beta-log(PF_L)/beta;
                    double betaDG_ii_L=beta*DG_ii_L;
                    sum=sum+P[i].nf_B_old*exp(-betaDG_ii_L);
                    sum_B=sum_B+P[i].nf_A_old*exp(-betaDG_ii_L);
                    
                    
                    P[i].nf_A_new=NA/(1+sum);
                    error=P[i].nf_A_old-P[i].nf_A_new;
                    
                    P[i].nf_B_new=NB/(1+sum_B);
                    error_B=P[i].nf_B_old-P[i].nf_B_new;
                
                    //cout << "old=" << P[i].nf_A_old << " new=" << P[i].nf_A_new << " error= " << fabs(error) << endl;
                    if (fabs(error) < 0.00000001 && fabs(error_B) < 0.00000001) Np_roots=Np_roots+1;
                    //if (fabs(error) < 0.0000001) Np_roots=Np_roots+1;
            
                    P[i].nf_A_old=P[i].nf_A_new;
                    P[i].nf_B_old=P[i].nf_B_new;
                }
                if (Np_roots==Np)   flag=false;
                Niteration=Niteration+1;
                //cout << "iteration= " << Niteration << endl;
            }
            //for (i=0; i<Np; i++) cout << "Free " << P[i].nf_A_old << " " << P[i].nf_B_old << endl;
            for (i=0; i<Np; i++)
            {
                double PF_L = 1/(P[i].omega_i*rho0);
                double DG_ii_L = DG0_L_star/beta-log(PF_L)/beta;
                //cout << DG_ii_L << endl;
                double betaDG_ii_L=beta*DG_ii_L;
                P[i].Pn[i].nl_d=P[i].nf_A_old*P[i].nf_B_old*exp(-betaDG_ii_L);
                
                //exit(0);
                for (j=0; j<P[i].maxnn; j++)
                {
                    index=P[i].Pn[j].nn_index;
                    
                    double PF = P[i].Pn[index].omega_ij/(P[i].omega_i*P[index].omega_i);
                    double DG_ij=DG0_star/beta-log(PF)/beta;
                    betaDG_ij=beta*DG_ij;
                    
                    P[i].Pn[index].nb_d=P[i].nf_A_old*P[index].nf_B_old*exp(-betaDG_ij)+
                        P[i].nf_B_old*P[index].nf_A_old*exp(-betaDG_ij);
                    cout << exp(-betaDG_ii_L) << " " << exp(-betaDG_ij) << " " <<  exp(-betaDG_ij)/exp(-betaDG_ii_L) << endl;
                }

            }
            
            double Free_energy_att=NA*log(P[0].nf_A_old/(double)NA) + NB*log(P[0].nf_B_old/(double)NB)
                +P[0].Pn[0].nl_d+(P[0].Pn[1].nb_d*z)/2.0;
            cout << P[0].nf_A_old+P[0].nf_B_old+P[0].Pn[0].nl_d*2+z*P[0].Pn[1].nb_d << endl;
                //exit(0);
            //for (i=0; i<Np; i++) cout << "bulk test " << P[i].nf_A_old << " " << P[i].nf_B_old << " " << P[i].Pn[i].nl_d << " " << P[0].Pn[1].nb_d*z << endl;
                
        for (i=0; i<Np; i++)
        {
           
            P[i].omega_i=0;
            P[i].Pn[i].nl_d=0;
            
            for (j=0; j<Np; j++)
            {
                P[i].Pn[j].nb_d=0;
               
                
                P[i].Pn[j].e_ij=0.0;
                P[i].Pn[j].e_ij_prime=0.0;
                P[i].Pn[j].omega_ij=0.0;
                P[i].Pn[j].omega_ij_prime=0.0;
                
            }   
        }    
                
            //exit(0);
            for (i=0; i<Np; i++)
            {   
                P[i].nf_A_old=rand()% NA+1;
                P[i].nf_B_old=rand()% NB+1;
                P[i].nf_A_new=P[i].nf_A_old;
                P[i].nf_B_new=P[i].nf_B_old;
                //cout << P[i].nf_A_old << endl;
            } 
            flag=true;    
            //cout << P[0].nf_A_old << " Random**** " << P[0].nf_B_old<< endl;
                
                /************************ Gas *****************/
                
            double gas_sum;
            double gas_sum_B;    
            while (flag==true && Niteration < 10000000)
            {
                Np_roots=0.0;
                for (i=0; i<Np-z; i++)
                {
                    gas_sum=0.0;
                    gas_sum_B=0.0;
                
                    double PF_L = 1/(omega_0*rho0);
                    double DG_ii_L = DG0_L_star/beta-log(PF_L)/beta;
                    double betaDG_ii_L=beta*DG_ii_L;
                    gas_sum=gas_sum+P[i].nf_B_old*exp(-betaDG_ii_L);
                    gas_sum_B=gas_sum_B+P[i].nf_A_old*exp(-betaDG_ii_L);
                    
                    
                    P[i].nf_A_new=NA/(1+gas_sum);
                    error=P[i].nf_A_old-P[i].nf_A_new;
                
                    P[i].nf_B_new=NB/(1+gas_sum_B);
                    error_B=P[i].nf_B_old-P[i].nf_B_new;
                
                    //cout << "old=" << P[i].nf_A_old << " new=" << P[i].nf_A_new << " error= " << fabs(error) << endl;
                    if (fabs(error) < 0.0000001 && fabs(error_B) < 0.0000001) Np_roots=Np_roots+1;
                    //if (fabs(error) < 0.0000001) Np_roots=Np_roots+1;
            
                    P[i].nf_A_old=P[i].nf_A_new;
                    P[i].nf_B_old=P[i].nf_B_new;
                }
            
            
                if (Np_roots==Np-z)   flag=false;
                Niteration=Niteration+1;
                //cout << "iteration= " << Niteration << endl;
            }
                //cout << P[0].nf_A_old << "**** " << P[0].nf_B_old<< endl; 
                //exit(0);
            for (i=0; i<Np-z; i++)
            {
                double PF_L = 1/(omega_0*rho0);
                double DG_ii_L = DG0_L_star/beta-log(PF_L)/beta;
                double betaDG_ii_L=beta*DG_ii_L;
                P[i].Pn[i].nl_d=P[i].nf_A_old*P[i].nf_B_old*exp(-betaDG_ii_L);
            } 
            
            
            double Gas_Free_energy_att=NA*log(P[0].nf_A_old/(double)NA) + NB*log(P[0].nf_B_old/(double)NB)
                +P[0].Pn[0].nl_d;
            double Free_energy_diff = Free_energy_att - Gas_Free_energy_att; 
                //cout << Free_energy_att << " " << Gas_Free_energy_att << " " << Free_energy_diff << endl;
            //if (N==90)myfile2 << DG0_star << " " << N << " " << Free_energy_diff << endl;
            //cout << "gas_test " << P[0].Pn[0].nl_d << " "  << Gas_Free_energy_att << endl;
           
            //cout << N << " " << Free_energy_diff << " " << target << endl; 
                 //exit(0);
            if (Free_energy_diff < target) {myfile2 << NA << " " << DG0_star << " " << 0 << endl; break;}

            //if (exp(Free_energy_diff-1)/CUB(delta/2.0) < rhog ) {myfile2 << rhog << " " << DG0_star << endl; break;}

            
            }
            //myfile2 << "\n"  << endl;
        }
    }
        
/****************************************************Force Calculations************************************************************************/
        
        
        
        
        //t = t + Delta_tB;
        cout << t << endl;
/**********************************************************************************************************************************************/   
        
        for (i=0; i<Np; i++)
        {
            
            P[i].omega_i=0;
            
            P[i].Pn[i].nl_d=0;
            
            for (j=0; j<Np; j++)
            {
                P[i].Pn[j].nb_d=0;
                
                P[i].Pn[j].e_ij=0.0;
                P[i].Pn[j].e_ij_prime=0.0;
                P[i].Pn[j].omega_ij=0.0;
                P[i].Pn[j].omega_ij_prime=0.0;
                
            }   
        }

        
/*********************************************************************************************************************************************/        
        
    }
    //myfile.close();
    //myfile1.close();
    myfile2.close();
    return 0;
}


void overlap(double r,double D,double D_,double *ov)
{
    //*ov=r*D*D_;
    //cout << "r= " << r << " D= " << D << " D_ " << D_ << endl;
    *ov=M_PI/(12.0*r) * SQ(D-r) * (SQ(r) + 2*r*D - 3*SQ(D_));
}
void overlap_prime(double r,double D,double D_,double *ov_prime)
{
    //*ov_prime= r*D*D_;
    *ov_prime  = (M_PI/(4.0*SQ(r))) * (D-r) * (D*SQ(D_) - D*SQ(r) + SQ(D_)*r - CUB(r));
    
}
void overlap_surface(double r,double D,double D_,double R2s, double *ov)
{
    *ov=M_PI/(3.0) * SQ(D-r) * (R2s + r + 2*D_);
}
void overlap_surface_prime(double r,double D,double D_,double R2s,double *ov_prime)
{
    *ov_prime=M_PI/(3.0) * (D-r) * (D-2*R2s-3*r-4*D_);
}
void normal(double *Rn)
{
    double x1, x2, w, y1, y2;
    do {
        //cout<< "FORCE CALCULATION2" << endl;
        x1 = 2.0 * rand()/(double)RAND_MAX - 1.0;
        x2 = 2.0 * rand()/(double)RAND_MAX - 1.0;
        w = x1 * x1 + x2 * x2;
        //cout<< "FORCE CALCULATION2=" << w << endl;
    } while ( w >= 1.0 );
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w; 
    *Rn=y1;
}

//void distance(Particles *P, double X, double Y, double Z,int i, double Lx, double Ly, double Lz,double *dx,double *dy,double *dz, double *d);
void distance_ads(double X, double Y, double Z,int i, double Lx, double Ly, double Lz,double *dx,double *dy,double *dz, double *d)
{
    //cout << P[i].x << " " << P[i].y << " " << P[i].z << endl;
    //exit(0);
    *dx = P[i].x-X;
    *dx = *dx - Lx * round(*dx/Lx);
            
    *dy = P[i].y-Y;
    *dy = *dy - Ly * round(*dy/Ly);
            
    *dz = P[i].z-Z;
    *dz = *dz - Lz * round(*dz/Lz);
            
    *d=sqrt( SQ(*dx) + SQ(*dy) + SQ(*dz));

}
void distance(int i, int j, double Lx, double Ly, double Lz,double *dx,double *dy,double *dz, double *d)
{
    //cout <<  "**" << P[i].omega_i << " " << P[i].Pn[j].omega_ij <<  endl;
    //cout << P[i].x << " " << P[i].y << " " << P[i].z << endl;
    //exit(0);
    *dx = P[i].x-P[j].x;
    *dx = *dx - Lx * round(*dx/Lx);
            
    *dy = P[i].y-P[j].y;
    *dy = *dy - Ly * round(*dy/Ly);
            
    *dz = P[i].z-P[j].z;
    *dz = *dz - Lz * round(*dz/Lz);
            
    *d=sqrt( SQ(*dx) + SQ(*dy) + SQ(*dz));

}
void overlapvolume(double r,double Rc,double L,double D1,double D1_,double D2,double D2_,double D12,double D12_,double *e_ij,double *w_ij,double *e_ij_prime,double *w_ij_prime)
{
    double ov, ov1, ov2, ov3, ov_prime, ov1_prime, ov2_prime, ov3_prime;
    if (r>(2*(Rc+L))) {
        *e_ij=0.0;                                       // corona-sphere overlap
        *w_ij=0.0;                                        //corona-corona overlap
        *e_ij_prime=0.0;
        *w_ij_prime=0.0;
    }
    else if (r<(2*(Rc+L)) && r>=2*Rc+1.9*L){
        *e_ij=0.0;
        *e_ij_prime=0.0;
        overlap(2*Rc+1.9*L,D1,D1_,&ov);
        overlap_prime(2*Rc+1.9*L,D1,D1_,&ov_prime);
        *w_ij=ov;
        *w_ij_prime=ov_prime;
    }
    else if (r<(2*Rc+1.9*L) && r>=2*Rc+L){
        *e_ij=0.0;
        *e_ij_prime=0.0;
                    
        overlap(r,D1,D1_,&ov);
        overlap_prime(r,D1,D1_,&ov_prime);
        *w_ij=ov;
        *w_ij_prime=ov_prime;
    }
    else if (r<2*Rc+L && r>2*Rc){
        overlap(r,D12,D12_,&ov);
        overlap_prime(r,D12,D12_,&ov_prime);
                    
        *e_ij=ov;
        *e_ij_prime=ov_prime;
                    
        overlap(r,D1,D1_,&ov);
        overlap_prime(r,D1,D1_,&ov_prime);
        *w_ij=ov-2.0**e_ij;
        *w_ij_prime=ov_prime-2.0**e_ij_prime;
    }
                
    else if (r<2*Rc){
        overlap(r,D12,D12_,&ov1);
        overlap(r,D2,D2_,&ov2);
        overlap_prime(r,D12,D12_,&ov1_prime);
        overlap_prime(r,D2,D2_,&ov2_prime);
        *e_ij=ov1-ov2;                                   //overlap(Rc+L,Rc)-overlap(Rc,Rc);
        *e_ij_prime=ov1_prime-ov2_prime;
                    
                    
        overlap(r,D1,D1_,&ov1);
        overlap(r,D12,D12_,&ov2);
        overlap(r,D2,D2_,&ov3);
        overlap_prime(r,D1,D1_,&ov1_prime);
        overlap_prime(r,D12,D12_,&ov2_prime);
        overlap_prime(r,D2,D2_,&ov3_prime);
        *w_ij=ov1 - 2*ov2 + ov3;                         //overlap(Rc+L,Rc+L)-2*overlap(Rc+L,Rc)+overlap(Rc,Rc);
        *w_ij_prime=ov1_prime - 2*ov2_prime + ov3_prime;            
    }
}
void overlapvolume_surface(double r,double Rc,double L,double R2s,double *e_is,double *e_si,double *w_is,double *e_is_prime,double *e_si_prime,double *w_is_prime)
{
    double ov, ov1, ov2, ov3, ov_prime, ov1_prime, ov2_prime, ov3_prime;
    double R1,R2,Dcs,dcs; 
    if (r>=Rc+L+R2s) {
        *e_is=0.0;                                       //corona-sphere overlap
        *w_is=0.0;                                       //corona-corona overlap
        *e_is_prime=0.0;
        *w_is_prime=0.0;
    }
    else if (r<(Rc+L+R2s) && r >= (Rc+R2s)){
        //cout << "***---" << r << " " << Rc+R2s << endl;
        R1=Rc+L; R2=R2s; Dcs=R1+R2; dcs=R1-R2;
        //cout << Dcs << " " << dcs << endl;
        *e_is=0.0;
        *e_is_prime=0.0;
        overlap_surface(r,Dcs,dcs,R2s,&ov);
        overlap_surface_prime(r,Dcs,dcs,R2s,&ov_prime);
        //cout << r << " " << Rc+R2s << " " << ov<< endl;
        *w_is=ov;
        *w_is_prime=ov_prime;
    }
    else if (r<Rc+R2s && r>=Rc){
        R1=Rc;R2=R2s;Dcs=R1+R2;dcs=R1-R2;
        overlap_surface(r,Dcs,dcs,R2s,&ov);
        overlap_surface_prime(r,Dcs,dcs,R2s,&ov_prime);
        *e_si=ov;       //surface linkes are blocked by the sphere
        *e_si_prime=ov_prime;   

        R1=Rc+L;R2=0.0;Dcs=R1+R2;dcs=R1-R2;
        overlap_surface(r,Dcs,dcs,R2,&ov);
        overlap_surface_prime(r,Dcs,dcs,R2,&ov_prime);
        *e_is=ov;       // colloids linkers are blocked by the surface
        *e_is_prime=ov_prime;
                
        R1=Rc+L;R2=R2s;Dcs=R1+R2;dcs=R1-R2; 
        overlap_surface(r,Dcs,dcs,R2s,&ov);
        overlap_surface_prime(r,Dcs,dcs,R2s,&ov_prime);
        *w_is=ov-*e_is-*e_si;
        *w_is_prime=ov_prime-*e_is_prime-*e_si_prime;
    }
}




