function P2 = P2_matrix(q)

    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    P2 = [    q3,    q2,    q1,    q0;
              q0,  - q1,    q2,  - q3;
            - q1,  - q0,    q3,    q2];
end