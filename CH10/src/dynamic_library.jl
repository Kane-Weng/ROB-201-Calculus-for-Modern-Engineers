function dyn_mod_planarSegway(q, dq)
    # DYN_MOD_PLANARSEGWAY
    # 2023-12-21 10:28:06
    #
    # Author: Grizzle
    #
    # Model NOTATION: D(q)ddq + C(q,dq)*dq + G(q) = B*tau 
    # The Robot Equations: From Lagrange's Equations of Motion
    #
    g, m_wh, m_pend, L, r_wh, J_pend, J_wh, J_rotor, N = modelParameters() 
    #
    # Variable names for the model
    theta, phi = q 
    dtheta, dphi = dq
    #
    D = zeros(2, 2)
    D[1, 1] = (J_pend*m_pend + J_pend*m_wh + (m_pend^2)*(r_wh^2) + (m_wh^2)*
                (r_wh^2) + 2.0J_wh*m_pend + 2.0J_wh*m_wh + (L^2)*(m_pend^2)*
                (cos(theta)^2) + (L^2)*(m_pend^2)*(sin(theta)^2) + 2.0m_pend*
                m_wh*(r_wh^2) + 2.0L*r_wh*(m_pend^2)*cos(theta) + 2.0L*m_pend*
                m_wh*r_wh*cos(theta)) / (m_pend + m_wh)
    D[1, 2] = 2.0J_wh + m_pend*(r_wh^2) + m_wh*(r_wh^2) + L*m_pend*r_wh*
                cos(theta)
    D[2, 1] = 2.0J_wh + m_pend*(r_wh^2) + m_wh*(r_wh^2) + L*m_pend*r_wh*
                cos(theta)
    D[2, 2] = 2.0J_wh + m_pend*(r_wh^2) + m_wh*(r_wh^2) + 2.0J_rotor*(N^2)
    #
    C = zeros(2, 2)
    C[1, 1] = (dtheta*(0.5(L^2)*(m_pend^2)*sin(2theta) - L*r_wh*(m_pend^2)*
                sin(theta) - (L^2)*(m_pend^2)*cos(theta)*sin(theta) - L*m_pend*
                m_wh*r_wh*sin(theta))) / (m_pend + m_wh)
    C[2, 1] = -L*dtheta*m_pend*r_wh*sin(theta)
    #
    G = zeros(2)
    G[1] = -L*g*m_pend*sin(theta)
    G[2] = 0
    #
    B = zeros(2, 1)
    B[2, 1] = 2N
    #
    JacG = zeros(2, 2)
    JacG[1, 1] = -L*g*m_pend*cos(theta)
    #
    JacComx = zeros(1,2)
    JacComx[1,1] = (m_pend*(r_wh + L*cos(theta)) + m_wh*r_wh) / (m_pend + m_wh)
    JacComx[1,2] = (m_pend*r_wh + m_wh*r_wh) / (m_pend + m_wh)
    #
    return (D=D, C=C, G=G, B=B, JacG=JacG, JacComx=JacComx)
end


# SEGWAY PARAMETERS
function modelParameters(outSelect=0)
    g = 9.81 # m/s^2
    m_wh = 6 # kg mass of 2 wheels 
    m_pend = 45+75 # kg mass of the person, platform, and the handlebar
    L = (45*.3 + 75*0.85)/(45+75) # m CoM for the person plus platform measured from
        # the axels of the wheels
    r_wh = 0.25 # m, Wheel radius
    J_pend = m_pend*L^2/6 # kg m^2 moment of inertia about the center of mass
    J_wh = 0.1 # kg m^2 moment of intertia of each wheel
    J_rotor = 1e-4 # kg m^2 per rotor, one for each wheel
    N = 50 # gear ratio between each motor rotor and wheel
    if outSelect > 0
        # named tuple
        return (g=g, m_wh=m_wh, m_pend=m_pend, L=L, r_wh=r_wh, J_pend=J_pend, 
            J_wh=J_wh, J_rotor=J_rotor, N=N)
    else # standard tuple
        return g, m_wh, m_pend, L, r_wh, J_pend, J_wh, J_rotor, N
    end
end
