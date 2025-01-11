#include "dynamicsolver.h"
#include <armadillo>
#include <iostream>

using namespace  std;
using namespace arma;

void DynamicSolver::SolveForwardDynamicTask()
{
    double g = 9.81;
    double ks = m_params->potForcesPar().ks;
    double xc = m_params->centGravity().Xc;
    double yc = m_params->centGravity().Yc;
    double a = m_params->geomPar().a;
    double b = m_params->geomPar().b;
    double c = m_params->geomPar().c;
    double R = m_params->geomPar().R;
    double L = m_params->geomPar().L;
    double Ix = m_params->inert().Ixx;
    double Iy = m_params->inert().Iyy;
    double m = m_params->potForcesPar().m;
    double zc = m_params->centGravity().Zc;
    double phi1 = m_jointsPos.phi1;
    double phi2 = m_jointsPos.phi2;
    double phi3 = m_jointsPos.phi3;
    double f0 = m_params->potForcesPar().f0;
    double C1a = L / (2 * b) * (cos(phi2) - cos(phi1));
    double C1b = L / (a + c)*(cos(phi3) - (cos(phi1) + cos(phi2))/2);
    double C2a = L / (4 * b)*(cos(phi1) * cos(phi1) - cos(phi2) * cos(phi2));
    double C2b = -L / (2 * (a + c))*(cos(phi3) * cos(phi3) - (cos(phi1) * cos(phi1) + cos(phi2) * cos(phi2))/2);
    double C2h = -L / (2 * (a + c))*(2 * a * sin(phi3)*sin(phi3) + c*(sin(phi1) * sin(phi1) + sin(phi2) * sin(phi2))/2);

    rowvec q = rowvec(3);
    q(0) = sin(phi1) - C1a;
    q(1) = sin(phi2) - C1a;
    q(2) = sin(phi3) - C1a;

    rowvec gamma = rowvec(3);
    gamma(0) = a / (2 * L) * C1b * C1b - C2b + b / L * C1a * C1b;
    gamma(1) = a / (2 * L) * C1b * C1b - C2b - b / L * C1a * C1b;
    gamma(2) = -c / (2 * L) * C1b * C1b - C2b;

    rowvec delta = rowvec(3);
    delta(0) = C2a - b / (2 * L) * C1a * C1a;
    delta(1) = C2a + b / (2 * L) * C1a * C1a;
    delta(2) = C2a;

    double kappa = R/L;
    double mu = (1/4) * C1a * (sin(phi2) - sin(phi1));
    double lambda = (1/4) * C1a * (2 * sin(phi3) - sin(phi1) - sin(phi2));

    rowvec e31 = rowvec(3);
    e31.ones();
    rowvec e31_ = {1, -1, 0};

    rowvec E_ = {a,a,-c};
    rowvec MuL = {-mu,mu,lambda};
    rowvec Sx0 = -kappa * C1b* e31 + pow(kappa,2)*gamma + pow(kappa,3)*(-mu*b/(a+c)*C1a + C1a*C2b*b/L*e31_ - C1b * C2h/L * e31);
    rowvec Sy0 = -kappa * q + pow(kappa,2)*delta + pow(kappa, 3)*(C1a*C2b*E_/L + C1a*C2h*e31/L + C1a*MuL);
    rowvec Sz0 = e31 - 1/2*pow(kappa, 2)*(C1b*C1b + q % q) + pow(kappa, 3) * (C1b * gamma + q % delta);

    rowvec F = rowvec(3);
    F(0) = 1 - 1/2*pow(kappa, 2)*pow(sin(phi1),2) - pow(kappa,3)*b/(2*L)*sin(phi1)*C1a*C1a;
    F(1) = 1 - 1/2*pow(kappa, 2)*pow(sin(phi2),2) + pow(kappa,3)*b/(2*L)*sin(phi2)*C1a*C1a;
    F(2) = 1 - 1/2*pow(kappa, 2)*pow(sin(phi3),2);

    rowvec D = rowvec(3);
    D = 1 - kappa*kappa*q + pow(kappa,3)*(C1b*gamma + q % delta);
    //D(0) = 1 - kappa^2/2*(q1^2+C1b^2) + kappa^3*(C1b*gamma(0) + q1*delta1);
    //D(1) = 1 - kappa^2/2*(q2^2+C1b^2) + kappa^3*(C1b*gamma2 + q2*delta2);
    //D(2) = 1 - kappa^2/2*(q3^2+C1b^2) + kappa^3*(C1b*gamma3 + q3*delta3);

    mat A0 = mat(3,3);
    A0(0,0) = F(0);
    A0(0,1) = F(1);
    A0(0,2) = F(2);
    A0(1,0) = b*D(0);
    A0(1,1) = -b*D(1);
    A0(1,2) = 0;
    A0(2,0) = a*D(0);
    A0(2,1) = a*D(1);
    A0(2,2) = -c*D(2);

    rowvec G = rowvec(3);
    G(0) = kappa * L * L * (sin(phi1) * (1 - kappa * cos(phi1)) - pow(kappa,2)/2*(b/L*C1a*C1a*cos(phi1) + pow(sin(phi1),3)) + pow(kappa,3)*b/(2*L)*C1a*C1a*(cos(2*phi1) + cos(phi1)*cos(phi2)));
    G(1) = kappa * L * L * (sin(phi2) * (1 - kappa * cos(phi2)) - pow(kappa,2)/2*(-b/L*C1a*C1a*cos(phi2) + pow(sin(phi2),3)) - pow(kappa,3)*b/(2*L)*C1a*C1a*(cos(2*phi2) + cos(phi1)*cos(phi2)));
    G(2) = kappa * L * L * (sin(phi3) * (1 - kappa * cos(phi3)) - pow(kappa,2)/2* pow(sin(phi3),3));

    rowvec zer3 = rowvec(3);
    zer3.zeros();
    mat Sxy0 = join_cols(Sy0, -Sx0);
    mat HorS0 = join_cols(zer3, Sxy0);
    mat A = A0 + HorS0*zc;
    colvec B = colvec(3);
    B(0) = (m * m_platformAccel.ddheight + ks * m_platformPos.height)/L + m * g * f0* L;
    B(1) = yc * m * g * L * (m_platformAccel.ddheight/g+1);
    B(2) = -xc * m * g * L *(m_platformAccel.ddheight/g+1);

    colvec k0 = colvec(3);
    solve(k0,A,B);
    colvec m0 = colvec(3);
    m0 = k0 % G.t();    // основной момент

    double mom_a2 = Ix * m_platformVel.droll * cos(m_platformPos.pitch);
    double mom_b2 = Iy * m_platformVel.dpitch;

    colvec B1 = {0, mom_a2/L, 0};
    colvec B2 = {0, 0, mom_b2/L};

    colvec X1 = colvec(3);
    colvec X2 = colvec(3);

    solve(X1, A, B1);
    solve(X2, A, B2);

    colvec md_a2 = colvec(3);
    colvec md_b2 = colvec(3);
    cout <<"k0 = " << k0.t() << endl;
    cout <<"m0 = " << m0.t() << endl;
    cout <<"Sx0 = " << Sx0 << endl;
    cout <<"Sy0 = " << Sy0 << endl;
    md_a2 = X1 % G.t();  // момент от второй производной угла крена
    md_b2 = X1 % G.t();  // момент от второй производной угла тангажа

    m_movingTorques.Q1 = m0(0) + md_a2(0) + md_b2(0);
    m_movingTorques.Q2 = m0(1) + md_a2(1) + md_b2(1);
    m_movingTorques.Q3 = m0(2) + md_a2(2) + md_b2(2);
}
