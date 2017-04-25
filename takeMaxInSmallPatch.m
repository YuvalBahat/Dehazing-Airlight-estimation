function maxInSmallPatch = takeMaxInSmallPatch(tLB,DCPpatchSize)
smallDCPpatchIndices = im2col(reshape(1:(sqrt(size(tLB,3))+DCPpatchSize-1)^2,ones(1,2)*sqrt(size(tLB,3))+DCPpatchSize-1),ones(1,2)*DCPpatchSize,'sliding');
maxInSmallPatch = permute(padarray(reshape(tLB,[size(tLB,1),ones(1,2)*sqrt(size(tLB,3))]),[0,ones(1,2)*floor(DCPpatchSize/2)],'replicate','both'),[2,3,1]);
maxInSmallPatch = max(maxInSmallPatch(bsxfun(@plus,smallDCPpatchIndices,shiftdim(0:size(maxInSmallPatch,3)-1,-1)*prod(partialSize(maxInSmallPatch,1:2)))));
maxInSmallPatch = shiftdim(maxInSmallPatch,2);
end