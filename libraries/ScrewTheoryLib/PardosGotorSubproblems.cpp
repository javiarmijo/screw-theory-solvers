// -*- mode:C++; tab-width:4; c-basic-offset:4; indent-tabs-mode:nil -*-

#include "ScrewTheoryIkSubproblems.hpp"

#include "ScrewTheoryTools.hpp"

#include <iostream>

using namespace roboticslab;

// -----------------------------------------------------------------------------

namespace
{
    KDL::Vector computeNormal(const MatrixExponential & exp1, const MatrixExponential & exp2)
    {
        KDL::Vector diff = exp2.getOrigin() - exp1.getOrigin();
        KDL::Vector normal = (exp1.getAxis() * diff) * exp1.getAxis();
        normal.Normalize();
        return vectorPow2(normal) * diff;
    }
}

// -----------------------------------------------------------------------------

PardosGotorOne::PardosGotorOne(const MatrixExponential & _exp, const KDL::Vector & _p)
    : exp(_exp),
      p(_p)
{}

// -----------------------------------------------------------------------------

bool PardosGotorOne::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector diff = k - f;
    double theta = KDL::dot(exp.getAxis(), diff);

    solutions = {{theta}};
    return true;
}

// -----------------------------------------------------------------------------

PardosGotorTwo::PardosGotorTwo(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const KDL::Vector & _p)
    : exp1(_exp1),
      exp2(_exp2),
      p(_p),
      crossPr2(exp2.getAxis() * exp1.getAxis()),
      crossPr2Norm(crossPr2.Norm())
{}

// -----------------------------------------------------------------------------

bool PardosGotorTwo::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector crossPr1 = exp2.getAxis() * (f - k);
    double crossPr1Norm = crossPr1.Norm();

    KDL::Vector c;

    if (KDL::dot(crossPr1, crossPr2) >= crossPr1Norm * crossPr2Norm)
    {
        c = k + (crossPr1Norm / crossPr2Norm) * exp1.getAxis();
    }
    else
    {
        c = k - (crossPr1Norm / crossPr2Norm) * exp1.getAxis();
    }

    double theta1 = KDL::dot(exp1.getAxis(), k - c);
    double theta2 = KDL::dot(exp2.getAxis(), c - f);

    solutions = {{theta1, theta2}};

    return true;
}

// -----------------------------------------------------------------------------

PardosGotorThree::PardosGotorThree(const MatrixExponential & _exp, const KDL::Vector & _p, const KDL::Vector & _k)
    : exp(_exp),
      p(_p),
      k(_k)
{}

// -----------------------------------------------------------------------------

bool PardosGotorThree::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector rhsAsVector = rhs * p - k;

    // std::cout<<"f = (" << f.x() << ", " << f.y() << ", " << f.z() << ")\n";
    // std::cout<<"pp = (" << (rhs * p).x() << ", " << (rhs * p).y() << ", " << (rhs * p).z() << ")\n";
    // std::cout<<"pg = (" << k.x() << ", " << k.y() << ", " << k.z() << ")\n";
    double delta = rhsAsVector.Norm();
    //std::cout << "delta = " << delta <<"\n";

    KDL::Vector diff = k - f;
    //std::cout<<"kmp = (" << diff.x() << ", " << diff.y() << ", " << diff.z() << ")\n";

    double dotPr = KDL::dot(exp.getAxis(), diff);
    //std::cout << "kmpp = " << dotPr << "\n";

    double sq2 = std::pow(dotPr, 2) - std::pow(diff.Norm(), 2) + std::pow(delta, 2);
    //std::cout << "root = " << sq2 << "\n";

    bool sq2_zero = KDL::Equal(sq2, 0.0);

    bool ret;
    //std::cout << "sq2 = " << sq2 <<"\n";
    if (!sq2_zero && sq2 > 0)
    {
        //std::cout <<"entra al if de pg3, por lo que ret es true\n";
        double sq = std::sqrt(std::abs(sq2));
        solutions = {{dotPr + sq}, {dotPr - sq}};
        ret = true;
    }
    else
    {
        //std::cout <<"entra al else de pg3\n";
        KDL::Vector proy = vectorPow2(exp.getAxis()) * diff;
        double norm = proy.Norm();
        solutions = {{norm}, {norm}};
        ret = sq2_zero;
        //if(ret) std::cout <<"pero ret es true\n";
    }

    //std::cout << "solutions_pg3 = {" << solutions[0][0] << ", " << solutions[1][0] << "}\n";

    return ret;
}

// -----------------------------------------------------------------------------

PardosGotorFour::PardosGotorFour(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const KDL::Vector & _p)
    : exp1(_exp1),
      exp2(_exp2),
      p(_p),
      n(computeNormal(exp1, exp2)),
      axisPow(vectorPow2(exp1.getAxis())) // same as exp2.getAxis()
{}

// -----------------------------------------------------------------------------

bool PardosGotorFour::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    // std::cout << "f = (" << f.x() << ", " << f.y() << ", " << f.z() << ")\n";
    // std::cout << "k = (" << k.x() << ", " << k.y() << ", " << k.z() << ")\n";
    // std::cout << "exp1.getOrigin() = (" << exp1.getOrigin().x() << ", " << exp1.getOrigin().y() << ", " << exp1.getOrigin().z() << ")\n";
    // std::cout << "exp2.getOrigin() = (" << exp2.getOrigin().x() << ", " << exp2.getOrigin().y() << ", " << exp2.getOrigin().z() << ")\n";
    // std::cout << "exp1.getAxis() = (" << exp1.getAxis().x() << ", " << exp1.getAxis().y() << ", " << exp1.getAxis().z() << ")\n";

    KDL::Vector u = f - exp2.getOrigin();
    KDL::Vector v = k - exp1.getOrigin();

    KDL::Vector u_p = u - axisPow * u;
    KDL::Vector v_p = v - axisPow * v;

    KDL::Vector c1 = exp1.getOrigin() + v - v_p;
    KDL::Vector c2 = exp2.getOrigin() + u - u_p;

    KDL::Vector c_diff = c2 - c1;
    // std::cout << "u = (" << u.x() << ", " << u.y() << ", " << u.z() << ")\n";
    // std::cout << "v = (" << v.x() << ", " << v.y() << ", " << v.z() << ")\n";
    // std::cout << "c1 = (" << c1.x() << ", " << c1.y() << ", " << c1.z() << ")\n";
    // std::cout << "c2 = (" << c2.x() << ", " << c2.y() << ", " << c2.z() << ")\n";
    // std::cout << "c_diff = (" << c_diff.x() << ", " << c_diff.y() << ", " << c_diff.z() << ")\n";
    // std::cout << "n = (" << n.x() << ", " << n.y() << ", " << n.z() << ")\n";
    bool samePlane = KDL::Equal(c_diff, n, 1e-4);

    if (!samePlane)
    {
        c_diff = n; // proyection of c_diff onto the perpendicular plane
        c1 = c2 - c_diff; // c1 on the intersecion of axis 1 and the normal plane to both axes
    }

    double c_norm = c_diff.Norm();
    double u_p_norm = u_p.Norm();
    double v_p_norm = v_p.Norm();

    double c_test = u_p_norm + v_p_norm - c_norm;
    bool c_zero = KDL::Equal(c_test, 0.0);

    // std::cout << "u_p_norm = " << u_p_norm << "\n";
    // std::cout << "v_p_norm = " << v_p_norm << "\n";
    // std::cout << "exp1.getAxis() = (" << exp1.getAxis().x() << ", " << exp1.getAxis().y() << ", " << exp1.getAxis().z() << ")\n";
    // std::cout << "exp2.getAxis() = (" << exp2.getAxis().x() << ", " << exp2.getAxis().y() << ", " << exp2.getAxis().z() << ")\n";

    if (!c_zero && c_test > 0.0 && u_p_norm > 0.0 && v_p_norm > 0.0)
    {
        //std::cout<<"if pg4\n";
        KDL::Vector omega_a = c_diff / c_norm;
        KDL::Vector omega_h = exp1.getAxis() * omega_a;

        double a = (std::pow(c_norm, 2) - std::pow(u_p_norm, 2) + std::pow(v_p_norm, 2)) / (2 * c_norm);
        double h = std::sqrt(std::abs(std::pow(v_p.Norm(), 2) - std::pow(a, 2)));

        KDL::Vector term1 = c1 + a * omega_a;
        KDL::Vector term2 = h * omega_h;

        KDL::Vector c = term1 + term2;
        KDL::Vector d = term1 - term2;

        KDL::Vector m1 = c - exp1.getOrigin();
        KDL::Vector m2 = c - exp2.getOrigin();

        KDL::Vector n1 = d - exp1.getOrigin();
        KDL::Vector n2 = d - exp2.getOrigin();

        KDL::Vector m1_p = m1 - axisPow * m1;
        KDL::Vector m2_p = m2 - axisPow * m2;

        KDL::Vector n1_p = n1 - axisPow * n1;
        KDL::Vector n2_p = n2 - axisPow * n2;

        double theta1_1 = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p), KDL::dot(m1_p, v_p));
        double theta2_1 = std::atan2(KDL::dot(exp2.getAxis(), u_p * m2_p), KDL::dot(u_p, m2_p));

        double theta1_2 = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p), KDL::dot(n1_p, v_p));
        double theta2_2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * n2_p), KDL::dot(u_p, n2_p));

        solutions = {
            {normalizeAngle(theta1_1), normalizeAngle(theta2_1)},
            {normalizeAngle(theta1_2), normalizeAngle(theta2_2)}
        };

        // if(samePlane) std::cout<<"samePlane\n";
        // std::cout << "c = (" << c.x() << ", " << c.y() << ", " << c.z() << ")\n";
        // std::cout << "d = (" << d.x() << ", " << d.y() << ", " << d.z() << ")\n";
        // std::cout << "c1 = (" << c1.x() << ", " << c1.y() << ", " << c1.z() << ")\n";
        // std::cout << "m1 = (" << m1.x() << ", " << m1.y() << ", " << m1.z() << ")\n";
        // std::cout << "m2 = (" << m2.x() << ", " << m2.y() << ", " << m2.z() << ")\n";
        // std::cout << "m1_p = (" << m1_p.x() << ", " << m1_p.y() << ", " << m1_p.z() << ")\n";
        // std::cout << "m2_p = (" << m2_p.x() << ", " << m2_p.y() << ", " << m2_p.z() << ")\n";
        // std::cout << "u = (" << u.x() << ", " << u.y() << ", " << u.z() << ")\n";
        // std::cout << "v = (" << v.x() << ", " << v.y() << ", " << v.z() << ")\n";
        // std::cout << "u_p = (" << u_p.x() << ", " << u_p.y() << ", " << u_p.z() << ")\n";
        // std::cout << "v_p = (" << v_p.x() << ", " << v_p.y() << ", " << v_p.z() << ")\n";
        // std::cout << "m1_p.Norm() = " << m1_p.Norm() << "\n";
        // std::cout << "m2_p.Norm() = " << m2_p.Norm() << "\n";
        // std::cout << "v_p_norm = " << v_p_norm << "\n";
        // std::cout << "u_p_norm = " << u_p_norm << "\n";

        return samePlane && KDL::Equal(m1_p.Norm(), v_p_norm);
    }
    else
    {
        //std::cout<<"else pg4\n";

        double theta1 = reference[0];
        double theta2 = reference[1];

        if (!KDL::Equal(v_p_norm, 0.0))
        {
            theta1 = std::atan2(KDL::dot(exp1.getAxis(), c_diff * v_p), KDL::dot(c_diff, v_p));
        }

        if (!KDL::Equal(u_p_norm, 0.0))
        {
            theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * c_diff), KDL::dot(-c_diff, u_p));
        }

        double normalized1 = normalizeAngle(theta1);
        double normalized2 = normalizeAngle(theta2);

        solutions = {
            {normalized1, normalized2},
            {normalized1, normalized2}
        };

        return samePlane && c_zero;
    }
}

// -----------------------------------------------------------------------------

PardosGotorFive::PardosGotorFive(const MatrixExponential & _exp, const MatrixExponential & _exp_next, const KDL::Vector & _p)
    : exp(_exp),
      exp_next(_exp_next),
      p(_p),
      axisPow(vectorPow2(exp.getAxis()))
{}

// -----------------------------------------------------------------------------

bool PardosGotorFive::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;
    KDL::Vector k_verify = rhs * p;
    // std::cout << "p = (" << p.x() << ", " << p.y() << ", " << p.z() << ")\n";
    // std::cout << "f = (" << f.x() << ", " << f.y() << ", " << f.z() << ")\n";
    // std::cout << "k = (" << k.x() << ", " << k.y() << ", " << k.z() << ")\n";
    // std::cout << "exp.getOrigin() = (" << exp.getOrigin().x() << ", " << exp.getOrigin().y() << ", " << exp.getOrigin().z() << ")\n";

    bool ret = true;

    //HACER QUE LA COORDENADA CORRESPONDIENTE AL EJE DE LA ROTACIÓN SEA LA MISMA PARA F(EQUIVALENTE A P) Y K
    for(int i=0; i < 3; i++)
    {
        if(exp.getAxis().data[i]!=0) k_verify.data[i]=f.data[i];
    }

    KDL::Vector u = f - exp.getOrigin();
    KDL::Vector v = k - exp.getOrigin();

    KDL::Vector u_w = axisPow * u;
    KDL::Vector v_w = axisPow * v;

    if (!(KDL::Equal(u_w, axisPow * (k_verify - exp.getOrigin())))) ret = false;

    KDL::Vector u_p = u - u_w;
    KDL::Vector v_p = v - v_w;

    double theta_k = reference[0];
    double theta_d = reference[0];

    if (!KDL::Equal(u_p.Norm(), 0.0) && !KDL::Equal(v_p.Norm(), 0.0))
    {
        // std::cout<<"calcula?\n";
        // std::cout << "v_p_norm = " << v_p.Norm() << "\n";
        // std::cout << "u_p_norm = " << u_p.Norm() << "\n";
        //theta_k = std::atan2(KDL::dot(exp.getAxis(), u_p * v_p), KDL::dot(u_p, v_p));
        double a = KDL::dot(exp.getAxis(), u_p * v_p);
        double b = KDL::dot(u_p, v_p);
        theta_k = std::atan2(a, b);
        theta_d= theta_k - KDL::PI;
    }
        //theta_d= theta_k - KDL::PI;

    //Ajuste PG5

    for(int i=0; i < 3; i++)
    {
        if(!(KDL::Equal(exp_next.getAxis().data[i],0)))
        {
            // std::cout<<"entra?\n";
            float x = dot(f - exp.getOrigin(), exp_next.getAxis());
            if(x != 0)
            {
                // std::cout<<"ajusta?\n";
                double d = f.data[i];//es la distancia que estará desplazado el plano con respecto al plano de movimiento

                //Recalcula los ángulso con el ajuste

                // std::cout <<"clamp(d / v_p.Norm()) = " << d / v_p.Norm() <<"\n";
                // std::cout <<"clamp(d / u_p.Norm()) = " << d / u_p.Norm() <<"\n";

                double sin1 = d / v_p.Norm();//acota el valor entre -1 y 1
                double sin2 = d / u_p.Norm();

                // std::cout << "theta_k = " << theta_k << "\n";
                // std::cout << "theta_d = " << theta_d << "\n";

                theta_k = theta_k - std::asin(sin1) + std::asin(sin2);
                theta_d = theta_d + std::asin(sin1) + std::asin(sin2);

                // theta_k = theta_k - d/v_p.Norm() + d/u_p.Norm();
                // theta_d = theta_k + d/v_p.Norm() + d/u_p.Norm();

                // std::cout << "theta_k = " << theta_k << "\n";
                // std::cout << "theta_d = " << theta_d << "\n";
            }
        }
    }

    //if (f==k) theta_d = theta_k; //SI P Y K SON IGAULES SOLO SACAMOS UNA SOLUCION, QUE SERA CERO. FUNCIONA ASI, PREGUNGTAR A BARTEK
    //theta_d = theta_k; //PARECE QUE NO SE DEBE CONSIDERAR LA SEGUNDA SOLUCION PARA OBTENER LAS SOLUCIONES ESPERADAS EN EL


    solutions = {{normalizeAngle(theta_k)}, {normalizeAngle(theta_d)}};

    // std::cout << "Solutions = {" << theta_k << "}, {" << theta_d << "}\n";
    // std::cout << "Solutions normalised = {" << normalizeAngle(theta_k) << "}, {" << normalizeAngle(theta_d) << "}\n";

    //return KDL::Equal(u_w, v_w);// && KDL::Equal(u_p.Norm(), v_p.Norm()); eso sería para pk1

    return ret;
}

// -----------------------------------------------------------------------------

PardosGotorSix::PardosGotorSix(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const KDL::Vector & _p)
    : exp1(_exp1),
      exp2(_exp2),
      p(_p),
      axesCross(exp1.getAxis() * exp2.getAxis()),
      axisPow1(vectorPow2(exp1.getAxis())),
      axisPow2(vectorPow2(exp2.getAxis())),
      axesDot(KDL::dot(exp1.getAxis(), exp2.getAxis()))
{}

// -----------------------------------------------------------------------------

bool PardosGotorSix::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector u = f - exp2.getOrigin();
    KDL::Vector v = k - exp1.getOrigin();

    KDL::Vector u_p = u - axisPow2 * u;
    KDL::Vector v_p = v - axisPow1 * v;

    KDL::Vector o2 = exp2.getOrigin() + axisPow2 * u;
    KDL::Vector o1 = exp1.getOrigin() + axisPow1 * v;

    double o2_dot = KDL::dot(exp2.getAxis(), o2);
    double o1_dot = KDL::dot(exp1.getAxis(), o1);

    KDL::Vector r3 = (exp1.getAxis() * (o1_dot - o2_dot * axesDot) + exp2.getAxis() * (o2_dot - o1_dot * axesDot)) / (1 - axesDot);

    MatrixExponential exp3(MatrixExponential::TRANSLATION, axesCross);
    PardosGotorThree pg3_2(exp3, r3, o2);
    PardosGotorThree pg3_1(exp3, r3, o1);

    Solutions pg3_1_sols, pg3_2_sols;

    bool pg3_1_ret = pg3_1.solve(KDL::Frame(v_p - (r3 - o1)), KDL::Frame::Identity(), pg3_1_sols);
    bool pg3_2_ret = pg3_2.solve(KDL::Frame(u_p - (r3 - o2)), KDL::Frame::Identity(), pg3_2_sols);

    KDL::Vector c2 = r3 + pg3_2_sols[0][0] * exp3.getAxis();
    KDL::Vector d2 = r3 + pg3_2_sols[1][0] * exp3.getAxis();

    KDL::Vector c1 = r3 + pg3_1_sols[0][0] * exp3.getAxis();
    KDL::Vector d1 = r3 + pg3_1_sols[1][0] * exp3.getAxis();

    bool ret = pg3_1_ret && pg3_2_ret;

    double theta1, theta2;

    if (c1 == c2)
    {
        KDL::Vector m2 = c2 - exp2.getOrigin();
        KDL::Vector m1 = c1 - exp1.getOrigin();

        KDL::Vector m2_p = m2 - axisPow2 * m2;
        KDL::Vector m1_p = m1 - axisPow1 * m1;

        theta1 = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p), KDL::dot(m1_p, v_p));
        theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * m2_p), KDL::dot(u_p, m2_p));
    }
    else if (d1 == d2)
    {
        KDL::Vector n2 = d2 - exp2.getOrigin();
        KDL::Vector n1 = d1 - exp1.getOrigin();

        KDL::Vector n2_p = n2 - axisPow2 * n2;
        KDL::Vector n1_p = n1 - axisPow1 * n1;

        theta1 = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p), KDL::dot(n1_p, v_p));
        theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * n2_p), KDL::dot(u_p, n2_p));
    }
    else
    {
        // these might be equal if `pg3_2_ret` is false
        KDL::Vector m1 = c2 - exp1.getOrigin();
        KDL::Vector n1 = d2 - exp1.getOrigin();

        if (m1.Norm() <= n1.Norm())
        {
            KDL::Vector m2 = c2 - exp2.getOrigin();

            KDL::Vector m2_p = m2 - axisPow2 * m2;
            KDL::Vector m1_p = m1 - axisPow1 * m1;

            theta1 = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p), KDL::dot(m1_p, v_p));
            theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * m2_p), KDL::dot(u_p, m2_p));
        }
        else
        {
            KDL::Vector n2 = d2 - exp2.getOrigin();

            KDL::Vector n2_p = n2 - axisPow2 * n2;
            KDL::Vector n1_p = n1 - axisPow1 * n1;

            theta1 = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p), KDL::dot(n1_p, v_p));
            theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * n2_p), KDL::dot(u_p, n2_p));
        }

        if(f != k) ret = false;
    }

    solutions = {{normalizeAngle(theta1), normalizeAngle(theta2)}};

    return ret;
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

PardosGotorSeven::PardosGotorSeven(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const MatrixExponential & _exp3, const KDL::Vector & _p)
    : exp1(_exp1),
      exp2(_exp2),
      exp3(_exp3),
      p(_p),
      n(computeNormal(exp1, exp2)),
      axisPow1(vectorPow2(exp1.getAxis())),
      axisPow2(vectorPow2(exp2.getAxis())),
      axesCross(exp1.getAxis() * exp2.getAxis()),
      axesCross_inverted(exp2.getAxis() * exp1.getAxis()),
      axesDot(KDL::dot(exp1.getAxis(), exp2.getAxis()))
{}

// -----------------------------------------------------------------------------

bool PardosGotorSeven::solve(const KDL::Frame &rhs, const KDL::Frame &pointTransform, const JointConfig &reference, Solutions &solutions) const {

    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;
    KDL::Vector u = f - exp3.getOrigin();
    KDL::Vector v = k - exp1.getOrigin();

    KDL::Vector u_p2 = u - axisPow2 * u;
    KDL::Vector v_p1 = v - axisPow1 * v;

    KDL::Vector o2 = exp2.getOrigin() + axisPow2 * (f - exp2.getOrigin());
    KDL::Vector o1 = exp1.getOrigin() + axisPow1 * v;

    double o2_dot = KDL::dot(exp2.getAxis(), o2);
    double o1_dot = KDL::dot(exp1.getAxis(), o1);

    KDL::Vector r4 = (exp1.getAxis() * (o1_dot - o2_dot * axesDot) + exp2.getAxis() * (o2_dot - o1_dot * axesDot)) / (1 - axesDot);
    KDL::Vector newAxes = KDL::dot(axesCross, o1 - o2) < 0.0 ? axesCross_inverted : axesCross;

    MatrixExponential exp4(MatrixExponential::TRANSLATION, newAxes);
    PardosGotorThree pg3_1(exp4, r4, o1);

    Solutions pg3_1_sols;
    if (!pg3_1.solve(KDL::Frame(v_p1 - (r4 - o1)), KDL::Frame::Identity(), pg3_1_sols)) return false;

    KDL::Vector c1 = r4 + (pg3_1_sols[0][0] * exp4.getAxis());
    KDL::Vector d1 = r4 + (pg3_1_sols[1][0] * exp4.getAxis());

    double theta_ck = reference[0], theta_dk = reference[1];
    PardosGotorFour pg4(exp2, exp3, f);

    Solutions pg4_c_sols, pg4_d_sols;
    bool pg4_ret_c = pg4.solve(KDL::Frame(c1 - f), KDL::Frame::Identity(), reference, pg4_c_sols);
    bool pg4_ret_d = pg4.solve(KDL::Frame(d1 - f), KDL::Frame::Identity(), reference, pg4_d_sols);

    if ((!pg4_ret_c && !pg4_ret_d) && (KDL::dot(axesCross, o1 - o2) == 0.0))
    {
        c1 = r4 + (-pg3_1_sols[0][0] * exp4.getAxis());
        d1 = r4 + (-pg3_1_sols[1][0] * exp4.getAxis());

        pg4_ret_c = pg4.solve(KDL::Frame(c1 - f), KDL::Frame::Identity(), reference, pg4_c_sols);
        pg4_ret_d = pg4.solve(KDL::Frame(d1 - f), KDL::Frame::Identity(), reference, pg4_d_sols);
    }

    if (pg4_ret_c && pg4_ret_d)
    {
        KDL::Vector m1 = c1 - exp1.getOrigin();
        KDL::Vector m1_p = m1 - axisPow1 * m1;

        KDL::Vector n1 = d1 - exp1.getOrigin();
        KDL::Vector n1_p = n1 - axisPow1 * n1;

        if (!KDL::Equal(v_p1.Norm(), 0.0))
        {
            theta_dk = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p1), KDL::dot(n1_p, v_p1));
            theta_ck = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p1), KDL::dot(m1_p, v_p1));
        }
    }
    else if (pg4_ret_c)
    {
        KDL::Vector m1 = c1 - exp1.getOrigin();
        KDL::Vector m1_p = m1 - axisPow1 * m1;

        if (!KDL::Equal(v_p1.Norm(), 0.0))
        {
            theta_ck = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p1), KDL::dot(m1_p, v_p1));
            theta_dk = theta_ck;
        }

        pg4_d_sols = pg4_c_sols;
    }
    else if (pg4_ret_d)
    {
        KDL::Vector n1 = d1 - exp1.getOrigin();
        KDL::Vector n1_p = n1 - axisPow1 * n1;

        if (!KDL::Equal(v_p1.Norm(), 0.0))
        {
            theta_dk = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p1), KDL::dot(n1_p, v_p1));
            theta_ck = theta_dk;
        }

        pg4_c_sols = pg4_d_sols;
    }
    else
    {
        return false;
    }

    solutions = {
        {theta_ck, pg4_c_sols[0][0], pg4_c_sols[0][1]},
        {theta_ck, pg4_c_sols[1][0], pg4_c_sols[1][1]},
        {theta_dk, pg4_d_sols[0][0], pg4_d_sols[0][1]},
        {theta_dk, pg4_d_sols[1][0], pg4_d_sols[1][1]}
    };

    return true;
}


// -----------------------------------------------------------------------------

PardosGotorEight::PardosGotorEight(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const MatrixExponential & _exp3, const KDL::Vector & _p, const int _firstID, const int _lastID,  const PoeExpression _poe)
    : exp1(_exp1),
      exp2(_exp2),
      exp3(_exp3),
      p(_p),
      firstID(_firstID),
      lastID(_lastID),
      poe(_poe),
      n(computeNormal(exp1, exp2)),
      axisPow(vectorPow2(exp1.getAxis())) // same as exp2.getAxis() and exp3
{}

// -----------------------------------------------------------------------------

bool PardosGotorEight::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions, const KDL::Frame & H_S_T, const KDL::JntArray & c_solutions, const KDL::Frame & H_S_T_0) const{
    
    KDL::Frame Hk;
    KDL::Frame frame_k;
    KDL::Frame Hp;
    KDL::Frame frame_p;

    // std::cout << "firstID= " << firstID << "\n";
    // std::cout << "lastID= " << lastID << "\n";

    for (int j = 0; j < firstID; j++)
    {
        const auto & exp = poe.exponentialAtJoint(j);
        frame_k = frame_k * exp.asFrame(c_solutions(j));
    }

    //for(int i = 0 ; i <= firstID; i++)
    //{
        Hk = frame_k.Inverse() * H_S_T;
    //}

    int x = lastID+1;

    for (int j = x; j < c_solutions.rows(); j++)
    {
        const auto & exp = poe.exponentialAtJoint(j);
        frame_p = frame_p * exp.asFrame(c_solutions(j));
    }

    //for(int i = 0 ; i <= firstID; i++)
    //{
        Hp = frame_p * H_S_T_0;
    //}


    //     std::cout << "\n Rotación Hk:" << std::endl; //H_S_T_0 ES HST0
    // for (int i = 0; i < 3; ++i) {
    //     std::cout << "  [ ";
    //     for (int j = 0; j < 3; ++j) {
    //         std::cout << Hk.M(i, j) << " ";
    //     }
    //     std::cout << "]" << std::endl;
    // }

    // std::cout << "Traslación Hk: [ "
    //           << Hk.p.x() << ", "
    //           << Hk.p.y() << ", "
    //           << Hk.p.z() << " ]\n" << std::endl;


    //     std::cout << "\n Rotación Hp:" << std::endl; //H_S_T_0 ES HST0
    // for (int i = 0; i < 3; ++i) {
    //     std::cout << "  [ ";
    //     for (int j = 0; j < 3; ++j) {
    //         std::cout << Hp.M(i, j) << " ";
    //     }
    //     std::cout << "]" << std::endl;
    // }

    // std::cout << "Traslación Hp: [ "
    //           << Hp.p.x() << ", "
    //           << Hp.p.y() << ", "
    //           << Hp.p.z() << " ]\n" << std::endl;


    // std::cout << "exp1.getOrigin() = ("<< exp1.getOrigin().x() << ", " << exp1.getOrigin().y() << ", " << exp1.getOrigin().z() << ") " << "\n"; 
    // std::cout << "exp2.getOrigin() = ("<< exp2.getOrigin().x() << ", " << exp2.getOrigin().y() << ", " << exp2.getOrigin().z() << ") " << "\n"; 
    // std::cout << "exp3.getOrigin() = ("<< exp3.getOrigin().x() << ", " << exp3.getOrigin().y() << ", " << exp3.getOrigin().z() << ") " << "\n"; 

    KDL::Vector f = Hp.p;
    KDL::Vector u3 = f - exp3.getOrigin();
    //std::cout << "u3 = ("<< u3.x() << ", " << u3.y() << ", " << u3.z() << ") " << "\n"; 
    KDL::Vector o3p = exp3.getOrigin() + axisPow *u3; 
   // std::cout << "o3p = ("<< o3p.x() << ", " << o3p.y() << ", " << o3p.z() << ") " << "\n";  
    KDL::Vector k_ = Hk.p;
    KDL::Vector o3k = Hk * Hp.Inverse() * o3p;
    //std::cout << "o3k = ("<< o3k.x() << ", " << o3k.y() << ", " << o3k.z() << ") " << "\n"; 
    
    Solutions pg4_sols;
    bool pg4_ret;

    //pg4

    KDL::Vector f_pg4 = o3p;
    KDL::Vector k_pg4 = o3k;

    KDL::Vector u = f_pg4 - exp2.getOrigin();
    KDL::Vector v = k_pg4 - exp1.getOrigin();

    KDL::Vector u_p = u - axisPow * u;
    KDL::Vector v_p = v - axisPow * v;

    // std::cout << "u_p = ("<< u_p.x() << ", " << u_p.y() << ", " << u_p.z() << ") " << "\n"; 
    // std::cout << "v_p = ("<< v_p.x() << ", " << v_p.y() << ", " << v_p.z() << ") " << "\n"; 

    KDL::Vector c1 = exp1.getOrigin() + v - v_p;
    KDL::Vector c2 = exp2.getOrigin() + u - u_p;

    // std::cout << "c1 = ("<< c1.x() << ", " << c1.y() << ", " << c1.z() << ") " << "\n"; 
    // std::cout << "c2 = ("<< c2.x() << ", " << c2.y() << ", " << c2.z() << ") " << "\n"; 

    KDL::Vector c_diff = c2 - c1;
    bool samePlane = KDL::Equal(c_diff, n, 1e-4);
    // if(samePlane) std::cout<<"sameplane\n";
    // if (!samePlane)
    {
        c_diff = n; // proyection of c_diff onto the perpendicular plane
        c1 = c2 - c_diff; // c1 on the intersecion of axis 1 and the normal plane to both axes
    }

    double c_norm = c_diff.Norm();
    double u_p_norm = u_p.Norm();
    double v_p_norm = v_p.Norm();

    double c_test = u_p_norm + v_p_norm - c_norm;
    // std::cout<< "c_test = " <<c_test << "\n";
    bool c_zero = KDL::Equal(c_test, 0.0);

    if (!c_zero && c_test > 0.0 && u_p_norm > 0.0 && v_p_norm > 0.0)
    {
        KDL::Vector omega_a = c_diff / c_norm;
        KDL::Vector omega_h = exp1.getAxis() * omega_a;

        double a = (std::pow(c_norm, 2) - std::pow(u_p_norm, 2) + std::pow(v_p_norm, 2)) / (2 * c_norm);
        double h = std::sqrt(std::abs(std::pow(v_p.Norm(), 2) - std::pow(a, 2)));

        KDL::Vector term1 = c1 + a * omega_a;
        KDL::Vector term2 = h * omega_h;

        KDL::Vector c = term1 + term2;
        KDL::Vector d = term1 - term2;

        KDL::Vector m1 = c - exp1.getOrigin();
        KDL::Vector m2 = c - exp2.getOrigin();

        KDL::Vector n1 = d - exp1.getOrigin();
        KDL::Vector n2 = d - exp2.getOrigin();

        KDL::Vector m1_p = m1 - axisPow * m1;
        KDL::Vector m2_p = m2 - axisPow * m2;

        KDL::Vector n1_p = n1 - axisPow * n1;
        KDL::Vector n2_p = n2 - axisPow * n2;

        double theta1_1 = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p), KDL::dot(m1_p, v_p));
        double theta2_1 = std::atan2(KDL::dot(exp2.getAxis(), u_p * m2_p), KDL::dot(u_p, m2_p));

        double theta1_2 = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p), KDL::dot(n1_p, v_p));
        double theta2_2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * n2_p), KDL::dot(u_p, n2_p));

        pg4_sols = {
            {normalizeAngle(theta1_1), normalizeAngle(theta2_1)},
            {normalizeAngle(theta1_2), normalizeAngle(theta2_2)}
        };
        
        // if(samePlane)
        // {
        //     std::cout << "m1_p.Norm() = "<< m1_p.Norm() << "\n";
        //     std::cout << "v_p_norm = "<< v_p_norm << "\n";
        // }

        pg4_ret = samePlane && KDL::Equal(m1_p.Norm(), v_p_norm);
    }
    else
    {
        double theta1 = reference[0];
        double theta2 = reference[1];

        if (!KDL::Equal(v_p_norm, 0.0))
        {
            theta1 = std::atan2(KDL::dot(exp1.getAxis(), c_diff * v_p), KDL::dot(c_diff, v_p));
        }

        if (!KDL::Equal(u_p_norm, 0.0))
        {
            theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * c_diff), KDL::dot(-c_diff, u_p));
        }

        double normalized1 = normalizeAngle(theta1);
        double normalized2 = normalizeAngle(theta2);

        pg4_sols = {
            {normalized1, normalized2},
            {normalized1, normalized2}
        };

        pg4_ret = samePlane && c_zero;
    }

    // std::cout <<"hola?\n";
    if(!pg4_ret) return false;
    // std::cout <<"hola?\n";

    //pk1

    KDL::Vector pk = o3k + (f - o3p);

    KDL::Vector u_w = axisPow * u;  
    KDL::Vector v_w = axisPow * v;

    KDL::Vector u_p_pk1 = pk - o3k;
    KDL::Vector v_p_pk1 = k_ - o3k;

    double theta123 = reference[0];

    if (!KDL::Equal(u_p_pk1.Norm(), 0.0) && !KDL::Equal(v_p_pk1.Norm(), 0.0))
    {
        theta123 = std::atan2(KDL::dot(exp3.getAxis(), u_p_pk1 * v_p_pk1), KDL::dot(u_p_pk1, v_p_pk1));
    }

    double theta3_1 = theta123 - pg4_sols[0][0] - pg4_sols[0][1];
    double theta3_2 = theta123 - pg4_sols[1][0] - pg4_sols[1][1];

    solutions = {
        {pg4_sols[0][0], pg4_sols[0][1], normalizeAngle(theta3_1)},
        {pg4_sols[1][0], pg4_sols[1][1], normalizeAngle(theta3_2)}
    };

    return KDL::Equal(u_w, v_w, 1e-2) && KDL::Equal(u_p_pk1.Norm(), v_p_pk1.Norm(), 1e-2);

} 

// -----------------------------------------------------------------------------

PardosGotorThreePadenKahanOne::PardosGotorThreePadenKahanOne(const MatrixExponential & _exp_pg3, const MatrixExponential & _exp_pk1, const KDL::Vector & _p, const KDL::Vector & _k)
    : exp_pg3(_exp_pg3),
      exp_pk1(_exp_pk1),
      p(_p),
      k(_k),
      axisPow(vectorPow2(exp_pk1.getAxis()))
{}

// -----------------------------------------------------------------------------

bool PardosGotorThreePadenKahanOne::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{   
    // std::cout << "exp_pg3.getAxis() = (" << exp_pg3.getAxis().x() << ", " << exp_pg3.getAxis().y() << ", " << exp_pg3.getAxis().z() << ")\n";
    // std::cout << "exp_pk1.getAxis() = (" << exp_pk1.getAxis().x() << ", " << exp_pk1.getAxis().y() << ", " << exp_pk1.getAxis().z() << ")\n";

    //  std::cout << "\nRotación RHS PG3_PK1:" << std::endl;
    // for (int i = 0; i < 3; ++i) {
    //     std::cout << "  [ ";
    //     for (int j = 0; j < 3; ++j) {
    //         std::cout << rhs.M(i, j) << " ";
    //     }
    //     std::cout << "]" << std::endl;
    // }

    // std::cout << "Traslación RHS PG3_PK1: [ "
    //           << rhs.p.x() << ", "
    //           << rhs.p.y() << ", "
    //           << rhs.p.z() << " ]\n" << std::endl;

    KDL::Vector f = pointTransform * p;
    KDL::Vector k2p = rhs * p;

    // std::cout << "k2p = (" << k2p.x() << ", " << k2p.y() << ", " << k2p.z() <<")\n";
    // std::cout << "p = (" << p.x() << ", " << p.y() << ", " << p.z() <<")\n";
    // std::cout << "f = (" << f.x() << ", " << f.y() << ", " << f.z() <<")\n";
    // std::cout << "k = (" << k.x() << ", " << k.y() << ", " << k.z() <<")\n";

    //PardosGotorThree pg3(exp_pg3, k2p, k);
    //PardosGotorThree pg3(exp_pg3, KDL::Vector(k2p.x(), k2p.y(), 0), KDL::Vector(k.x(), k.y(), 0));//EN MATLAB PARDOS LO HACE ASI
    Solutions pg3_sols;

    //bool pg3_ret = pg3.solve(rhs, pointTransform, reference, pg3_sols);
    //bool pg3_ret = pg3.solve(KDL::Frame::Identity(), pointTransform, reference, pg3_sols);

    // std::cout <<"holaaaaaaaaaaaaaaa? \n";

    //if(!pg3_ret) return false;

    //PG3

    //KDL::Vector f = pointTransform * p;
    KDL::Vector rhsAsVector = f - k;

    // std::cout<<"f = (" << f.x() << ", " << f.y() << ", " << f.z() << ")\n";
    // std::cout<<"pp = (" << (rhs * p).x() << ", " << (rhs * p).y() << ", " << (rhs * p).z() << ")\n";
    // std::cout<<"pg = (" << k.x() << ", " << k.y() << ", " << k.z() << ")\n";
    // std::cout<<"rhsAsVector = (" << rhsAsVector.x() << ", " << rhsAsVector.y() << ", " << rhsAsVector.z() << ")\n";
    double delta = rhsAsVector.Norm();
    // std::cout << "delta = " << delta <<"\n";

    
    //KDL::Vector diff = k - k2p;
    KDL::Vector diff = KDL::Vector(k.x(), k.y(), 0.0) - KDL::Vector(k2p.x(), k2p.y(), 0.0);
    // std::cout<<"kmp = (" << diff.x() << ", " << diff.y() << ", " << diff.z() << ")\n";

    double dotPr = KDL::dot(exp_pg3.getAxis(), diff);
    // std::cout << "kmpp = " << dotPr << "\n";

    double sq2 = std::pow(dotPr, 2) - std::pow(diff.Norm(), 2) + std::pow(delta, 2);
    // std::cout << "root = " << sq2 << "\n";

    bool sq2_zero = KDL::Equal(sq2, 0.0);

    bool pg3_ret;
    // std::cout << "sq2 = " << sq2 <<"\n";
    if (!sq2_zero && sq2 > 0)
    {
        // std::cout <<"entra al if de pg3, por lo que ret es true\n";
        double sq = std::sqrt(std::abs(sq2));
        pg3_sols = {{normalizeAngle(dotPr + sq)}, {normalizeAngle(dotPr - sq)}};
        pg3_ret = true;
    }
    else
    {
        // std::cout <<"entra al else de pg3\n";
        KDL::Vector proy = vectorPow2(exp_pg3.getAxis()) * diff;
        double norm = proy.Norm();
        pg3_sols = {{normalizeAngle(norm)}, {normalizeAngle(norm)}};
        pg3_ret = sq2_zero;
        // if(pg3_ret) std::cout <<"pero ret es true\n";
    }

    // std::cout << "solutions_pg3 = {" << pg3_sols[0][0] << ", " << pg3_sols[1][0] << "}\n";

    if(!pg3_ret) return false;


    //FIN PG3



    // std::cout<<"llega aqui??  \n";

    KDL::Vector k2 = k2p + exp_pg3.getAxis() * pg3_sols[1][0];
    // std::cout << "k2 = (" << k2.x() << ", " << k2.y() << ", " << k2.z() <<")\n";
    
    // std::cout << "pg3_sols[0][0] = " << pg3_sols[0][0] << "\n"; 
    // std::cout << "pg3_sols[0][1] = " << pg3_sols[0][1] << "\n"; 
    // std::cout << "pg3_sols[1][0] = " << pg3_sols[1][0] << "\n"; 
    // std::cout << "pg3_sols[1][1] = " << pg3_sols[1][1] << "\n"; 

    // std::cout << "k2p = (" << k2p.x() << ", " << k2p.y() << ", " << k2p.z() <<")\n";
    // std::cout << "k2 = (" << k2.x() << ", " << k2.y() << ", " << k2.z() <<")\n";

    // std::cout << "exp_pk1.getOrigin() = (" << exp_pk1.getOrigin().x() << ", " << exp_pk1.getOrigin().y() << ", " << exp_pk1.getOrigin().z() << ")\n";

    KDL::Vector u = f - exp_pk1.getOrigin();
    KDL::Vector v = k2 - exp_pk1.getOrigin();

    // std::cout << "u = (" << u.x() << ", " << u.y() << ", " << u.z() <<")\n";
    // std::cout << "v = (" << v.x() << ", " << v.y() << ", " << v.z() <<")\n";

    /*
    for(int i=0; i < 3; i++)
    {
        if(exp_pk1.getAxis().data[i]!=0) v.data[i]=u.data[i];
    }
    */

    //v.data[1] = u.data[1];

    // std::cout << "u = (" << u.x() << ", " << u.y() << ", " << u.z() <<")\n";
    // std::cout << "v = (" << v.x() << ", " << v.y() << ", " << v.z() <<")\n";

    // std::cout << "v.data[1] = " << v.data[1] << "\n";
    // std::cout << "u.data[1] = " << u.data[1] << "\n";

    KDL::Vector u_w = axisPow * u;
    KDL::Vector v_w = axisPow * v;

    KDL::Vector u_p = u - u_w;
    KDL::Vector v_p = v - v_w;

    double theta = reference[0];

    // std::cout << "u_p.Norm() = " << u_p.Norm() <<"\n";
    // std::cout << "v_p.Norm() = " << v_p.Norm() <<"\n";
    // std::cout << "u_w = (" << u_w.x() << ", " << u_w.y() << ", " << u_w.z() <<")\n";
    // std::cout << "v_w = (" << v_w.x() << ", " << v_w.y() << ", " << v_w.z() <<")\n";

    if (!KDL::Equal(u_p.Norm(), 0.0) && !KDL::Equal(v_p.Norm(), 0.0))
    {
        theta = std::atan2(KDL::dot(exp_pk1.getAxis(), u_p * v_p), KDL::dot(u_p, v_p));
    }

    theta = normalizeAngle(theta);

    // std::cout<<"tetha = " << theta << "\n";
    // std::cout<<"tetha normalised = " << normalizeAngle(theta) << "\n";
    // std::cout<<"--tetha = " << -theta << "\n";
    // std::cout<<"-tetha normalised = " << normalizeAngle(-theta) << "\n";

    //theta=reference[0];//PUEDE SER UTIL

/*
    if(KDL::Equal(abs(theta), abs(reference[0]))) solutions = {{theta}, {-theta}};
    else solutions = {{theta}, {theta}};
*/

    solutions = {{theta}, {-theta}};

    /*
    double theta2= theta - KDL::PI;
    solutions = {{theta}, {theta2}};
    */
   return KDL::Equal(u_p.Norm(), v_p.Norm());//KDL::Equal(u_w, v_w) && COMENTADO

}

// -----------------------------------------------------------------------------

Algebraic_UR::Algebraic_UR(int _j1, int _j2) //RECIBE ENTEROS CON LAS POSICIONES DE LAS ARTICULACIONES RESUELTAS
    : j1(_j1),
      j2(_j2)
{}

// -----------------------------------------------------------------------------

bool Algebraic_UR::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions, const KDL::Frame & H_S_T_0, const KDL::JntArray & c_solutions, const KDL::Frame & quitar) const
{
    //TIENE QUE RECIBIR SOLUTIONS Y H_S_T_0

    //     std::cout << "\nRotación RHS ALGEBRAIC:" << std::endl;
    // for (int i = 0; i < 3; ++i) {
    //     std::cout << "  [ ";
    //     for (int j = 0; j < 3; ++j) {
    //         std::cout << rhs.M(i, j) << " ";
    //     }
    //     std::cout << "]" << std::endl;
    // }

    // std::cout << "Traslación RHS ALGEBRAIC: [ "
    //           << rhs.p.x() << ", "
    //           << rhs.p.y() << ", "
    //           << rhs.p.z() << " ]\n" << std::endl;

    double nx = H_S_T_0.M(0, 0);
    double ny = H_S_T_0.M(1, 0);
    double ox = H_S_T_0.M(0, 1);
    double oy = H_S_T_0.M(1, 1);


    // std::cout <<"nx = " << nx << "\n";
    // std::cout <<"ny = " << ny << "\n";
    // std::cout <<"ox = " << ox << "\n";
    // std::cout <<"oy = " << oy << "\n";
    // std::cout << "c_solutions(j1)= " << c_solutions(j1) <<"\n";
    // std::cout << "c_solutions(j2)= " << c_solutions(j2) <<"\n";
    // std::cout << "s1= " << sin(c_solutions(j1)) <<"\n";
    //  std::cout << "s1_rounded= " << std::round(std::sin(c_solutions(j1)) * 1e4) / 1e4 <<"\n";
    // std::cout << "c1= " << cos(c_solutions(j1)) <<"\n";
    // std::cout << "s5= " << sin(c_solutions(j2)) <<"\n";

    double theta = reference[0];

    //if(sin(c_solutions(j2)) != 0) theta = atan2((ox * sin(c_solutions(j1)) - oy * cos(c_solutions(j1))) / sin(c_solutions(j2)) ,(ny * cos(c_solutions(j1)) - nx * sin(c_solutions(j1)) / sin(c_solutions(j2))));
    
    /*
    if(sin(c_solutions(j2)) != 0) 
    {
        theta = std::round(
            std::atan2(
                std::round((
                    std::round(ox * 1e4) / 1e4 * std::round(std::sin(c_solutions(j1)) * 1e4) / 1e4 -
                    std::round(oy * 1e4) / 1e4 * std::round(std::cos(c_solutions(j1)) * 1e4) / 1e4
                ) / (std::round(std::sin(c_solutions(j2)) * 1e4) / 1e4) * 1e4) / 1e4,

                std::round((
                    std::round(ny * 1e4) / 1e4 * std::round(std::cos(c_solutions(j1)) * 1e4) / 1e4 -
                    std::round(nx * 1e4) / 1e4 * std::round(std::sin(c_solutions(j1)) * 1e4) / 1e4
                ) / (std::round(std::sin(c_solutions(j2)) * 1e4) / 1e4) * 1e4) / 1e4
            ) * 1e4
        ) / 1e4;
    }
    */

    // std::cout<<"sin(c_solutions(j2)) = " << sin(c_solutions(j2)) << "\n";

    if(sin(c_solutions(j2)) != 0)
    {
        //constexpr double epsilon = 1e-10;

        double s1 = std::sin(c_solutions(j1));
        double c1 = std::cos(c_solutions(j1));
        double s2 = std::sin(c_solutions(j2)); 

        double A = ox * s1 - oy * c1;
        double B = ny * c1 - nx * s1;

        double C = A / s2;
        double D = B / s2;

        // std::cout <<"s1 = " << s1 << "\n";
        // std::cout <<"c1 = " << c1 << "\n";
        // std::cout <<"s2 = " << s2 << "\n";
        // std::cout <<"A = " << A << "\n";
        // std::cout <<"B = " << B << "\n";
        // std::cout <<"C = " << C << "\n";
        // std::cout <<"D = " << D << "\n";

        theta = std::atan2(C, D); 
    }
    

    /*
    std::cout << "solutions = ( ";
    for(int i=0; i < c_solutions.rows(); i ++)
    {
        std::cout << c_solutions(i)<< "  ";
    }
    std::cout << "\n";
    */
    // std::cout << "tetha= " << theta <<"\n";
    // std::cout << "tetha normalised= " << normalizeAngle(theta) <<"\n";

    //theta=reference[0];//PUEDE SER UTIL

    solutions = {{normalizeAngle(theta)}};

    return true; //???????????
}

// -----------------------------------------------------------------------------
