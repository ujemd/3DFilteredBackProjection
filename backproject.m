function f = backproject(q, rVec, phiVec, xVec, yVec, intpol)

[xGrid, yGrid] = meshgrid(xVec, yVec);
f = zeros(size(xGrid));
Nr = length(rVec);
Nphi = length(phiVec);
deltaR = rVec(2) - rVec(1);

switch intpol
  case 'nearest'
    for phiIx = 1:Nphi
      phi = phiVec(phiIx);
      proj = q(:, phiIx);
      r = xGrid * cos(phi) + yGrid * sin(phi);
      rIx1 = round(r / deltaR + (Nr + 1) / 2);
      f = f + proj(rIx1);
    end
  case 'linear'
    for phiIx = 1:Nphi
      phi = phiVec(phiIx);
      proj = q(:, phiIx);
      r = xGrid * cos(phi) + yGrid * sin(phi);
      rIx = r / deltaR + (Nr + 1) / 2;
      rIx1 = min(ceil(rIx),Nr);
      rIx2 = min(floor(rIx),Nr);
      w1 = rIx1-rIx;
      w2 = rIx-rIx2;
      f = f + w2.*proj(rIx1) + w1.*proj(rIx2);
    end
  otherwise
    error('Unknown interpolation.');
end

f = f * pi / (Nphi*deltaR);