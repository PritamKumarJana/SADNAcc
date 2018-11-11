# include <stdio.h>
# include <sys/types.h>
# include <math.h>
# include <stdlib.h>
# include <time.h>
# include <fstream>
# include <iostream>
# define maxNp 100
# define maxN 50
# define maxh 30

# define SQ(x) ((x)*(x))
# define CUB(x) ((x)*(x)*(x))
using namespace std;

void overlapvolume(double,double,double,double,double,double,double,double,double ,double *,double *,double *,double *);
void overlapvolume_surface(double,double,double,double ,double *,double *,double *,double *,double *,double *);
void overlap(double,double,double,double *);
void overlap_prime(double,double,double,double *);
void overlap_surface(double,double,double,double,double *);
void overlap_surface_prime(double,double,double,double,double *);

typedef struct 
{
	double nf_A,nf_B;			//Maximum number of free linkers, Maximum number of neighbor
    double nf_A_new,nf_B_new;
    double nbS_AC;
    double nbS_d;
}Particles;
Particles P[maxNp];

int main()
{
    int N;
    int NA,NB,NsC;
    double beta=1.0;
    int hmax_=30;
    double F[maxh]={0.0};
    int hmax;
    //double Lx=(2.0*11.0)/pow(3,1.0/4.0);
	//double Ly=(2.0*11.0)/pow(3,1.0/4.0);
    double Lx=11.0;
	double Ly=11.0;
    //double A=0.5*Lx*sqrt(3)/2.0*Lx;
    //double Nreceptor=100.0;
    //double receptor_density=(Nreceptor/A);
    cout << 0.5*Lx*sqrt(3)/2.0*Lx<< endl;
    //exit(0);
	double Lz=75.0;
    double nfS_C;
    double nfS_C_new;
    double Rc=5.0;
	double L=1.0;
    //double L=2.0;
    //double L=0.5
    double R2s=1.0;
     
    bool flag;
    double rho0=1.0;
    
    //int seed =time(NULL)*getpid();
    //int seed = 14;
    //srand (seed);
    //cout << "seed=" << seed << endl;
    
    
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
    double omega_0=4*M_PI*SQ(Rc)*L;
    double e_ij_sum,omega_i;
    
    double e_is;
    double e_is_prime;
    double e_si;
    double e_si_prime;
    double w_is;
    double w_is_prime;
    double e_si_sum;
    double omega_S;
    double omega_S_0=Lx*Ly*R2s;
    e_si_sum=0.0;
    int h;
    double omega_ij,omega_is;
    double error=0.000001;
    double PF,PF_L,PF_surface;
    double DG_ij,betaDG_ij,qb_colloid,DG_ij_surface,DG0_star,DG0_star_surface,betaDG_ij_surface,qb_surface,DG_ii_L,DG0_star_L, betaDG_ii_L,ql;
    double nf_A_bulk,nf_B_bulk,nf_A_bulk_new,nf_B_bulk_new;
    int counter_h;
    //double rhog=0.01;
    double vol=(4.0/3.0)*M_PI*CUB(5.0);
    //double rhog=0.1/vol;
    double rho_ideal=4.1e-06;
    ofstream outfile;
    //ofstream outfile1;
    outfile.open ("SMFig4c.txt");
    //outfile1.open ("Phase_diagram_surface_fcc-1.dump");
    int DG0max=15;
    int Nmax=1000;
    int rmax=1;  
    double e=2.71828;
    for (int counter_r=0; counter_r<rmax; counter_r++)
    {
        double r=11.0;
        overlapvolume(r,Rc,L,D1,D1_,D2,D2_,D12,D12_,&e_ij,&w_ij,&e_ij_prime,&w_ij_prime);
        r = 5.7;
        overlapvolume_surface(r,Rc,L,R2s,&e_is,&e_si,&w_is,&e_is_prime,&e_si_prime,&w_is_prime);
        omega_ij=w_ij;
        omega_is=w_is;
        for (int counter_DG0=-5; counter_DG0 < DG0max; counter_DG0++)
        {
            DG0_star=-9.0; //+ counter_DG0; 
            DG0_star_L=DG0_star; 
            DG0_star_surface=DG0_star+counter_DG0; //-1.0;// -10.0;
            //rho_ideal=0.000000000001 + 0.00000001*counter_DG0;
            for (int counter_N=0; counter_N<Nmax; counter_N=counter_N+1)
            {
                //counter_N=DG0_star-DG0_star_surface;
                // DG0_star_surface=DG0_star+counter_N;
                //N=2+counter_N;
                N=80;
                NA=N/2;
                NB=N/2;
                NsC=215; //+counter_N*10.0;
                rho_ideal=0.000000001 + 0.00000001*counter_N;

                nf_A_bulk = rand()% NA+1;
                nf_B_bulk = rand()% NB+1;
                
                cout << nf_A_bulk << " " << nf_B_bulk << endl;
                PF_L = 1/(omega_0*rho0);
                
                DG_ii_L = DG0_star_L/beta-log(PF_L)/beta;
                
                betaDG_ii_L=beta*DG_ii_L;
                ql=exp(-betaDG_ii_L);
                flag=true;
                while (flag==true)
                {            
                    nf_A_bulk_new=NA/(1 + nf_B_bulk*ql);
                    nf_B_bulk_new=NB/(1 + nf_A_bulk*ql);
                    double diff_A=fabs(nf_A_bulk_new-nf_A_bulk);
                    double diff_B=fabs(nf_B_bulk_new-nf_B_bulk);
                    
                    if (diff_A<error && diff_B<error) flag=false;
                    nf_A_bulk=nf_A_bulk_new;
                    nf_B_bulk=nf_B_bulk_new;
                }
               
                //double betaFreeEnergy_bulk=NA*log((nf_A_bulk*nf_B_bulk)/(NA*NB))+(NA-nf_A_bulk);
                double FreeEnergy_bulk=(1/beta)*(NA*log((nf_A_bulk*nf_B_bulk)/(NA*NB))+(NA-nf_A_bulk));
                
                //cout << nf_A_bulk << " " << nf_B_bulk << " " << FreeEnergy_bulk << endl;
                //cout << nf_A_bulk << " " << nf_B_bulk << " " << (NA-nf_A_bulk) << endl;
                
                double Grand_partition=0.0;
                for (hmax=0; hmax<=hmax_; hmax++)   //hmax_
                {
                    for (h=0; h<hmax; h++)
                    {
                        
                        P[h].nf_A=rand()% NA+1;
                        P[h].nf_B=rand()% NB+1;
                    }
//                    cout << "test1" << endl;
//                    exit(0);
                    nfS_C=rand()% NsC+1;
                    flag=true;
                    while (flag==true)
                    {
                        if (hmax==0) flag=false;
                        counter_h=0;
                        for (h=0; h<hmax; h++)
                        {
                            if (h==0)
                            {
                                if (hmax==1) e_ij_sum=6*e_ij+e_is;
                                else e_ij_sum=9*e_ij+e_is;
                                omega_i=omega_0-e_ij_sum;
                                
                                e_si_sum=e_si;
                                omega_S=omega_S_0-e_si_sum;
                                
                                PF = omega_ij/(omega_i*omega_i);
                                
                                DG_ij=DG0_star/beta-log(PF)/beta;
                                betaDG_ij=beta*DG_ij;
                                qb_colloid=exp(-betaDG_ij);
                                
                                
                                PF_surface = omega_is/(omega_i*omega_S);
                                //cout << log(PF) << " "<< log(PF_surface) << endl; 
                                
                                DG_ij_surface=DG0_star_surface/beta-log(PF_surface)/beta;
                                betaDG_ij_surface = beta*DG_ij_surface;
                                qb_surface = exp(-betaDG_ij_surface);
                                
                                //qb_surface=5.0;
                                
                                PF_L = 1/(omega_i*rho0);
                                DG_ii_L = DG0_star_L/beta-log(PF_L)/beta;
                                betaDG_ii_L=beta*DG_ii_L;
                                ql=exp(-betaDG_ii_L);
                                
                                //ql=20.0;
                                
                                nfS_C_new=NsC/(1+P[h].nf_A*qb_surface);
                                if (hmax==1)
                                {
                                    P[h].nf_A_new=NA/(1 + 6*P[h].nf_B*qb_colloid + P[h].nf_B*ql + nfS_C*qb_surface);
                                    P[h].nf_B_new=NB/(1 + 6*P[h].nf_A*qb_colloid + P[h].nf_A*ql);
                                }
                                else {
                                    P[h].nf_A_new=
                                        NA/(1 + 6*P[h].nf_B*qb_colloid + P[h].nf_B*ql + nfS_C*qb_surface + 3*P[h+1].nf_B*qb_colloid);
                                    P[h].nf_B_new=
                                        NB/(1 + 6*P[h].nf_A*qb_colloid + P[h].nf_A*ql + 3*P[h+1].nf_A*qb_colloid);
                                }
                                
                            }
                            else if (h>0 && h!=(hmax-1))
                            {
                                e_ij_sum=12*e_ij;
                                omega_i=omega_0-e_ij_sum;
                            
                                PF = omega_ij/(omega_i*omega_i);
                                DG_ij=DG0_star/beta-log(PF)/beta;
                                betaDG_ij=beta*DG_ij;
                                qb_colloid=exp(-betaDG_ij);
                                //qb_colloid=5.0;
                            
                                PF_L = 1/(omega_i*rho0);
                                DG_ii_L = DG0_star_L/beta-log(PF_L)/beta;
                                betaDG_ii_L=beta*DG_ii_L;
                                ql=exp(-betaDG_ii_L);
                                //ql=20.0;
                                
                                P[h].nf_A_new=
                                    NA/(1 + 6*P[h].nf_B*qb_colloid + P[h].nf_B*ql + 3*P[h-1].nf_B*qb_colloid + 3*P[h+1].nf_B*qb_colloid);
                                P[h].nf_B_new=
                                    NB/(1 + 6*P[h].nf_A*qb_colloid + P[h].nf_A*ql + 3*P[h-1].nf_A*qb_colloid + 3*P[h+1].nf_A*qb_colloid);
                            }
                            else if (h==hmax-1)
                            {
                                e_ij_sum=9*e_ij;
                                omega_i=omega_0-e_ij_sum;
                            
                                PF = omega_ij/(omega_i*omega_i);
                                DG_ij=DG0_star/beta-log(PF)/beta;
                                betaDG_ij=beta*DG_ij;
                                qb_colloid=exp(-betaDG_ij);
                                //qb_colloid=5.0;
                                
                                PF_L = 1/(omega_i*rho0);
                                DG_ii_L = DG0_star_L/beta-log(PF_L)/beta;
                                betaDG_ii_L=beta*DG_ii_L;
                                ql=exp(-betaDG_ii_L);
                                
                                //ql=20.0;
                                
                                P[h].nf_A_new=
                                    NA/(1 + 6*P[h].nf_B*qb_colloid + P[h].nf_B*ql + 3*P[h-1].nf_B*qb_colloid);
                                P[h].nf_B_new=
                                    NB/(1 + 6*P[h].nf_A*qb_colloid + P[h].nf_A*ql + 3*P[h-1].nf_A*qb_colloid);
                            }
                            //cout << ql << " " << qb_colloid << " " << qb_surface << " " << qb_colloid/ql << endl; 
                            double diff_A=fabs(P[h].nf_A_new-P[h].nf_A);
                            double diff_B=fabs(P[h].nf_B_new-P[h].nf_B);
                            double diff_C=fabs(nfS_C_new-nfS_C);
                            if (diff_A<error && diff_B<error && diff_C<error) {counter_h=counter_h+1;}
                        
                            P[h].nf_A=P[h].nf_A_new;
                            P[h].nf_B=P[h].nf_B_new;
                            nfS_C=nfS_C_new;
                        
                            if (counter_h==hmax) flag=false;  
                        }
                    }
                    //cout << "****" << P[0].nf_B*P[0].nf_A*qb_colloid << " " << P[0].nf_A*nfS_C*qb_surface << " " << P[0].nf_A*P[0].nf_B*ql << endl;
                    
                    double sum=0.0;
                    for (h=0; h<hmax; h++) 
                    {
                        sum=sum+NA*log((P[h].nf_A*P[h].nf_B)/(NA*NB)) + (2*NA-P[h].nf_A-P[h].nf_B)/2.0;
                    }
                    //cout << "\n" << endl;
                    if (hmax==0) F[hmax]=1.0;
                    else 
                    {
                        double FreeEnergy=(1/beta)*(NsC*log(nfS_C/(double)NsC) + (NsC-nfS_C)/2.0 + sum);
                        double Free_energy_diff = FreeEnergy-hmax*FreeEnergy_bulk; 
                        double delta=1.0;
                        F[hmax]=exp(-beta*Free_energy_diff)*pow((e*rho_ideal*CUB(delta/2)),hmax);
                    }
                    
                    
                    //outfile1 << hmax << " " << NsC*log(nfS_C/(double)NsC) + (NsC-nfS_C)/2.0 << " " << nfS_C/(double)NsC << endl;
                    
                    //F[hmax]=exp(-betaFreeEnergy+hmax*betaFreeEnergy_bulk);//*pow((rhog*CUB(L/2)),hmax);
                    
                    Grand_partition=Grand_partition+F[hmax];
                    
                    
//                    sum_betaFreeEnergy=sum_betaFreeEnergy+betaFreeEnergy;
//                    cout << betaFreeEnergy << endl;
//                    outfile << hmax << " " << exp(-betaFreeEnergy+hmax*betaFreeEnergy_bulk) << endl;
                }
                double gamma=0.0;
                for (hmax=0; hmax<=hmax_; hmax++) 
                {
                    //outfile << N << " " << hmax << " " << F[hmax]/Grand_partition << endl;
                    //if (hmax==0) gamma=F[hmax]/Grand_partition;
                    //else 
                    gamma=gamma+hmax*F[hmax]/Grand_partition;
                    //outfile << hmax << " " << F[hmax]/(double)hmax << endl;
                    F[hmax]=0.0;
                }
                //outfile << DG0_star << " " << N << " " << gamma << endl;
                //outfile << DG0_star << " " << gamma << " " << DG0_star_surface << endl;
                //outfile << DG0_star << " " << N << " " << gamma << " " << DG0_star_surface << endl;
                //outfile << rho_ideal << " " << gamma << " " << N << " " << DG0_star_surface << endl;
                outfile << DG0_star << " " << DG0_star_surface << "  " << rho_ideal*e << " " << gamma << " " << DG0_star-DG0_star_surface << endl;
                //outfile << (rho_ideal) << " " << N << " " << gamma << endl;
                //outfile << "\n" << endl;
            }
            outfile << "\n" ;
        }
    }
    outfile.close();
}
void overlap(double r,double D,double D_,double *ov)
{
    *ov=M_PI/(12.0*r) * SQ(D-r) * (SQ(r) + 2*r*D - 3*SQ(D_));
}

void overlap_prime(double r,double D,double D_,double *ov_prime)
{
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
    else if (r<2*Rc+L && r>=2*Rc){
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
    else if (r<(Rc+L+R2s) && r >= (Rc+L+0.9*R2s)){
        //cout << "***---" << r << " " << Rc+R2s << endl;
        R1=Rc+L; R2=R2s; Dcs=R1+R2; dcs=R1-R2;
        //cout << Dcs << " " << dcs << endl;
        r=Rc+L+0.9*R2s;
        *e_is=0.0;
        *e_is_prime=0.0;
        overlap_surface(r,Dcs,dcs,R2,&ov);
        overlap_surface_prime(r,Dcs,dcs,R2s,&ov_prime);
        //cout << r << " " << Rc+R2s << " " << ov<< endl;
        *w_is=ov;
        *w_is_prime=ov_prime;
    }
    else if (r < (Rc+L+0.9*R2s) && r >= (Rc+L)){
        //cout << "***---" << r << " " << Rc+R2s << endl;
        R1=Rc+L; R2=R2s; Dcs=R1+R2; dcs=R1-R2;
        //cout << Dcs << " " << dcs << endl;
        *e_is=0.0;
        *e_is_prime=0.0;
        overlap_surface(r,Dcs,dcs,R2,&ov);
        overlap_surface_prime(r,Dcs,dcs,R2s,&ov_prime);
        //cout << r << " " << Rc+R2s << " " << ov<< endl;
        *w_is=ov;
        *w_is_prime=ov_prime;
    }
    
    else if (r<Rc+R2s && r>=Rc){
        R1=Rc;R2=R2s;Dcs=R1+R2;dcs=R1-R2;
        overlap_surface(r,Dcs,dcs,R2,&ov);
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
//    else if (){
//        R1=Rc+0.75*L;R2=0.0;Dcs=R1+R2;dcs=R1-R2;
//        overlap_surface(r,Dcs,dcs,R2,&ov);
//        overlap_surface_prime(r,Dcs,dcs,R2,&ov_prime);
//        *e_is=ov;       // colloids linkers are blocked by the surface
//        *e_is_prime=ov_prime;
//        
//        R1=Rc;R2=R2s*0.75;Dcs=R1+R2;dcs=R1-R2;
//        overlap_surface(r,Dcs,dcs,R2,&ov);
//        overlap_surface_prime(r,Dcs,dcs,R2,&ov_prime);
//        *e_si=ov;       //surface linkes are blocked by the sphere
//        *e_si_prime=ov_prime;
//        
//    }
}

