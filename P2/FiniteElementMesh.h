#pragma once

#include "AnimatedTetrahedronMesh.h"

#include "CGVectorWrapper.h"
#include "CGSystemWrapper.h"
#include "ConjugateGradient.h"

#include <Eigen/Dense>
#include <string> //debugging purposes

template<class T>
struct FiniteElementMesh : public AnimatedTetrahedronMesh<T>
{
    using Base = AnimatedTetrahedronMesh<T>;
    using Base::m_meshElements;
    using Base::m_particleX;
    using Vector3 = typename Base::Vector3;
    using Matrix33 = Eigen::Matrix< T , 3 , 3>;

    int m_nFrames;
    int m_subSteps;
    T m_frameDt;
    T m_stepDt;
    T m_stepEndTime;

    const T m_density;
    const T m_mu;
    const T m_lambda;
    const T m_rayleighCoefficient;

    std::vector<T> m_particleMass;
    std::vector<Vector3> m_particleV;
    std::vector<Matrix33> m_DmInverse;
    std::vector<T> m_restVolume;
    
    FiniteElementMesh(const T density, const T mu, const T lambda, const T rayleighCoefficient)
        :m_density(density), m_mu(mu), m_lambda(lambda), m_rayleighCoefficient(rayleighCoefficient)
    {}

    void initializeUndeformedConfiguration()
    {
        // Initialize rest shape and particle mass (based on constant density)
        m_particleMass.resize(m_particleX.size(), T()); // Initialize all particle masses to zero
        for(const auto& element: m_meshElements)
        {
            Matrix33 Dm;
            for(int j = 0; j < 3; j++)
                Dm.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            T restVolume = .5 * Dm.determinant();
            if(restVolume < 0)
                throw std::logic_error("Inverted element" + std::to_string(restVolume));
            m_DmInverse.emplace_back(Dm.inverse());
            m_restVolume.push_back(restVolume);
            T elementMass = m_density * restVolume;
            for(const int v: element)
                m_particleMass[v] += (1./3.) * elementMass;
        }
    }
    
    void addElasticForce(std::vector<Vector3>& f) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Linear Elasticity
            Matrix33 Ds;
            for(int j = 0; j < 2; j++)
                Ds.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            Matrix33 F = Ds * m_DmInverse[e];

            Matrix33 strain = .5 * (F + F.transpose()) - Matrix33::Identity();
            Matrix33 P = 2. * m_mu * strain + m_lambda * strain.trace() * Matrix33::Identity();

            Matrix33 H = -m_restVolume[e] * P * m_DmInverse[e].transpose();
            
            for(int j = 0; j < 2; j++){
                f[element[j+1]] += H.col(j);
                f[element[0]] -= H.col(j);
            }
        }
    }

    void addProductWithStiffnessMatrix(std::vector<Vector3>& w, std::vector<Vector3>& f, const T scale) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Linear Damping
            Matrix33 Ds_dot;
            for(int j = 0; j < 2; j++)
                Ds_dot.col(j) = w[element[j+1]]-w[element[0]];
            Matrix33 F_dot = Ds_dot * m_DmInverse[e];

            Matrix33 strain_rate = .5 * (F_dot + F_dot.transpose());
            Matrix33 P_damping = scale * (2. * m_mu * strain_rate + m_lambda * strain_rate.trace() * Matrix33::Identity());

            Matrix33 H_damping = m_restVolume[e] * P_damping * m_DmInverse[e].transpose();
            
            for(int j = 0; j < 2; j++){
                f[element[j+1]] += H_damping.col(j);
                f[element[0]] -= H_damping.col(j);
            }
        }
    }

    void simulateSubstep()
    {
        using FEMType = FiniteElementMesh<T>;        

        const int nParticles = m_particleX.size();

        // Save last time step velocities

        std::vector<Vector3> lastV = m_particleV;
        
        // Construct initial guess for next-timestep
        //   Velocities -> Same as last timestep
        //   Positions -> Using Forward Euler
        
        for(int p = 0; p < nParticles; p++)
            m_particleX[p] += m_stepDt * m_particleV[p];

        // Overwrite boundary conditions with desired values

        setBoundaryConditions();
        
        // Solve for everything else using Conjugate Gradients

        std::vector<Vector3> dx(nParticles, Vector3::Zero());
        std::vector<Vector3> rhs(nParticles, Vector3::Zero());
        std::vector<Vector3> q(nParticles, Vector3::Zero());
        std::vector<Vector3> s(nParticles, Vector3::Zero());
        std::vector<Vector3> r(nParticles, Vector3::Zero());
        CGVectorWrapper<Vector3> dxWrapper(dx);
        CGVectorWrapper<Vector3> rhsWrapper(rhs);
        CGVectorWrapper<Vector3> qWrapper(q);
        CGVectorWrapper<Vector3> sWrapper(s);
        CGVectorWrapper<Vector3> rWrapper(r);
        CGSystemWrapper<Vector3, FEMType> systemWrapper(*this);
        
        addElasticForce(rhs);
        for(int p = 0; p < nParticles; p++)
            rhs[p] += (m_particleMass[p] / m_stepDt) * (lastV[p] - m_particleV[p]);
        addProductWithStiffnessMatrix(m_particleV, rhs, -m_rayleighCoefficient);
        clearConstrainedParticles(rhs);

        ConjugateGradient<T>::Solve(systemWrapper,
            dxWrapper, rhsWrapper, qWrapper, sWrapper, rWrapper,
            1e-4, 50);

        // Apply corrections to positions and velocities

        const T oneOverDt = T(1.) / m_stepDt;
        for(int p = 0; p < nParticles; p++){
            m_particleX[p] += dx[p];
            m_particleV[p] += oneOverDt * dx[p];
        }
    }

    void simulateFrame(const int frame)
    {
        m_stepDt = m_frameDt / (T) m_subSteps;

        for(int step = 1; step <= m_subSteps; step++){
            m_stepEndTime = m_frameDt * (T) (frame-1) + m_stepDt * (T) step;
            simulateSubstep();
        }
    }

    virtual void clearConstrainedParticles(std::vector<Vector3>& x) {}
    virtual void setBoundaryConditions() {}
};

