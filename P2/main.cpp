//#include "AnimatedTetrahedronMesh.h" //borrowed from Demo3D
#include "FiniteElementMesh.h"

#include <Eigen/Dense>

#include <map>
#include <cstdlib>

template<class T>
struct LatticeMesh : public FiniteElementMesh<T>
{
    using Base = FiniteElementMesh<T>;

    // from AnimatedTetrahedonMesh
    using Base::m_meshElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;
    using Vector3 = typename Base::Vector3;
 
    // from FiniteElementsMesh
    using Base::m_particleV;
    using Base::initializeUndeformedConfiguration;
    using Base::m_stepEndTime;


    std::array<int, 3> m_cellSize; // dimensions in grid cells
    int m_boxSize; // box size in grid cells
    T m_gridDX;

    std::vector<Vector3> m_particleUndeformedX;
    std::vector<int> m_leftHandleIndicies;
    std::vector<int> m_rightHandleIndicies;
    //Add vectors for the deformation
    Vector3 m_leftHandleVelocity;
    Vector3 m_rightHandleVelocity;   
    //Vector3 m_topHandleVelocity;


    std::vector<std::array<int, 3>> m_activeCells; // Marks the "active" cells in the lattice
    std::map<std::array<int, 3>, int> m_activeNodes; // Maps the "active" nodes to their particle index

    //MN add constructor
   LatticeMesh():Base(1.e2, 1., 4., .05)
	{
		m_leftHandleVelocity = Vector3(0, 0., 0.);  //just one edge moves
		m_rightHandleVelocity = Vector3(.2, 0., 0.);
		//m_topHandleVelocity = Vector3(.2, 0., 0.);
  	 }

    void initialize()
    {
        initializeUSD("NewbyP2.usda");

        // Activate cells in a box in the center of our space

        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
        for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++)
        for(int cell_k = 0; cell_k < m_cellSize[2]; cell_k++){  

           // if(r <= m_radius * m_radius)
           if(abs(cell_i - m_cellSize[0]/2) <= m_boxSize/2)
		    if(abs(cell_j - m_cellSize[0]/2) <= m_boxSize/2)
			     if(abs(cell_k - m_cellSize[0]/2) <= m_boxSize/2)
	   				m_activeCells.push_back(std::array<int, 3>{cell_i, cell_j, cell_k});

        }

        std::cout << "Created a model including " << m_activeCells.size() << " lattice cells" <<std::endl;

        // Create (uniquely numbered) particles at the node corners of active cells

        for(const auto& cell: m_activeCells){
            std::array<int, 3> node;
            for(node[0] = cell[0]; node[0] <= cell[0]+1; node[0]++)
            for(node[1] = cell[1]; node[1] <= cell[1]+1; node[1]++)
            for(node[2] = cell[2]; node[2] <= cell[2]+1; node[2]++){
                auto search = m_activeNodes.find(node);
                if(search == m_activeNodes.end()){ // Particle not yet created at this lattice node location -> make one
                    m_activeNodes.insert({node, m_particleX.size()});
                    m_particleX.emplace_back(m_gridDX * T(node[0]), m_gridDX * T(node[1]), m_gridDX * T(node[2]));
                }
            }
        }
        std::cout << "Model contains " << m_particleX.size() << " particles" << std::endl;

        // Make tetrahedra out of all active cells (6 tetrahedra per cell)

        for(const auto& cell: m_activeCells){
            int vertexIndices[2][2][2];
            for(int i = 0; i <= 1; i++)
            for(int j = 0; j <= 1; j++)
            for(int k = 0; k <= 1; k++){
                std::array<int, 3> node{cell[0] + i, cell[1] + j, cell[2] + k};
                auto search = m_activeNodes.find(node);
                if(search != m_activeNodes.end())
                    vertexIndices[i][j][k] = search->second;
                else
                    throw std::logic_error("particle at cell vertex not found");
            }

            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][0], vertexIndices[1][1][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][1], vertexIndices[1][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][1], vertexIndices[1][1][1], vertexIndices[0][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][1], vertexIndices[0][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][0], vertexIndices[0][1][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][0], vertexIndices[0][1][0], vertexIndices[1][1][1]});
        }
        
        // Perform the USD-specific initialization of topology & particles
        // (this will also create a boundary *surface* to visualuze

        initializeTopology();
        initializeParticles();

	// Also resize the velocities to match
        m_particleV.resize(m_particleX.size(), Vector3::Zero());

        // Initialize rest shape matrices and particle mass
        initializeUndeformedConfiguration();

        // Also record rest shape
        m_particleUndeformedX = m_particleX;

	//Find edges of the box. We assume it's a cube here. 
	//this would need to be more complext to handle rectangles (one for each dim) 
	float minBox = m_gridDX * (m_cellSize[0]/2. - m_boxSize/2.) ;
	float maxBox = minBox + m_gridDX * m_boxSize;
//TODO Remove debugging junk
std::cout << "minBox: " << minBox << "  MaxBox: " << maxBox << std::endl;
	for(const auto& particle: m_particleX){
		//left handle is the edge with x/y min
                if(particle[0] == minBox)
		if(particle[1] == minBox){
			//std::cout << "Particle: " << particle[0] << particle[1]<< particle[2] << std::endl;
			m_leftHandleIndicies.push_back(gridToParticleID(particle[0], particle[1], particle[2]));
		}
		//Right handle is the edge with x/y min
                if(particle[0] >= maxBox)
                if(particle[1] >= maxBox){
                        //std::cout << "Max Particle: " << particle[0] << particle[1]<< particle[2] << std::endl;
                        m_rightHandleIndicies.push_back(gridToParticleID(particle[0], particle[1], particle[2]));
                }
	}
	    //std::cout << "left handle contains " << m_leftHandleIndicies.size() << " particles" << std::endl;
	    //std::cout << "Right handle contains " << m_rightHandleIndicies.size() << " particles" << std::endl;
	    
	// Check particle indexing in mesh

        for(const auto& element: m_meshElements)
            for(const auto vertex: element)
                if(vertex < 0 || vertex >= m_particleX.size())
                    throw std::logic_error("mismatch between mesh vertex and particle array");
    }

    void initializeDeformation()
    {
        // No need to apply any deformation; this example is driven by moving handles
    }

    void clearConstrainedParticles(std::vector<Vector3>& x) override
    {
     //   for(const auto v: m_leftHandleIndicies)
     //       x[v] = Vector3::Zero();
    //    for(const auto v: m_rightHandleIndicies)
    //        x[v] = Vector2::Zero();
    }

    void setBoundaryConditions() override
    {
        T effectiveTime = std::min<T>(m_stepEndTime, 1.0);
	//TODO: Add top if desired.
        for(const auto v: m_leftHandleIndicies){
            m_particleX[v] = m_particleUndeformedX[v]; // + effectiveTime * m_leftHandleVelocity;
            m_particleV[v] = m_leftHandleVelocity;
        }
       // for(const auto v: m_rightHandleIndicies){
       //     m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_rightHandleVelocity;
        //    m_particleV[v] = m_rightHandleVelocity;
       // }
    }

    
private:
    inline int gridToParticleID(const int i, const int j, const int k) const { return i * (m_cellSize[1]+1) + j * (m_cellSize[2]+1) + k; }

};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 20, 20, 20 };
    simulationMesh.m_boxSize = 8;
    simulationMesh.m_gridDX = 0.1;

    simulationMesh.m_nFrames = 300;
    simulationMesh.m_subSteps = 1;
    simulationMesh.m_frameDt = 0.1;

    // Initialize the simulation example
    simulationMesh.initialize();

    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);
    
    //iterate through frames of animation
    for(int frame=1; frame <= simulationMesh.m_nFrames; frame++){
	simulationMesh.simulateFrame(frame);
	simulationMesh.writeFrame(frame);
    }


    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}

