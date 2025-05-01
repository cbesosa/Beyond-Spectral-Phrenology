function pdfISI = calculate_pdf_ISI(spkT)
                
    ISI = spkT(2:end)-spkT(1:end-1);
    recursiveISI = nan(length(spkT)*length(ISI)/2,1);
    pointer = 1;
    for ik = 1:length(ISI)
        recursiveISI(pointer:pointer+length(spkT)-ik-1) = spkT(1+ik:end)-spkT(1:end-ik);
        pointer = pointer +length(spkT)-ik;
    end
    [numCounts, binEdges] = histcounts(1./ISI, 0:1:512);
    pdfISI.Recip = numCounts./length(ISI);  % reciprocal.
    pdfISI.binCenters = (binEdges(1:end-1) + binEdges(2:end))/2;

    [numCounts, binEdges] = histcounts(1./recursiveISI, 0:1:512);
    pdfISI.ReRecip = numCounts./length(recursiveISI);  % reciprocal of recursive ISI.
    
    % nan problem:
    pdfISI.Recip(isnan(pdfISI.Recip)) = 0;
    pdfISI.ReRecip(isnan(pdfISI.ReRecip)) = 0;
end
    