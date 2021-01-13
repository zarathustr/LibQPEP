function J = J_pnp_loss(image_pt, world_pt, K, R, t)
numPoints = size(image_pt, 1);
J = 0;
for i = 1 : numPoints
    world_point = [world_pt(i, :), 1]; % homogeneous coordinates
    image_point = image_pt(i, :);

    cameraMatrix = [R; t.'] * K;
    projectedPoint = world_point * cameraMatrix;
    
    if(~isnumeric(R))
        image_point = image_point .* projectedPoint(3);
        projectedPoint = projectedPoint(1 : 2);
    else
        projectedPoint = projectedPoint(1 : 2) / projectedPoint(3);
    end
    d = image_point - projectedPoint;
    J = J + d * d';
end
J = J ./ numPoints;
end