function [chains]=getChains(ijlist,chains)

max_len = 10;

if nargin<2,
    chains = cell(0);
    chainz = cell(0);
    
    upre = unique(ijlist(:,1));
    for i=1:length(upre)
        tmp = getChains(ijlist,{upre(i)});
        for j=1:length(tmp)
            chains{(length(chains)+1)}=tmp{j};
        end
    end
else
    chain_ext = cell(0);
    upos = unique(ijlist(ijlist(:,1)==chains{1}(end),2));
    for i=1:length(upos)
        if ~any(chains{1}==upos(i)) & length(chains{1})<max_len
            tmp = getChains(ijlist,{[chains{1} upos(i)]});
            for j=1:length(tmp)
                chains{(length(chains)+1)}=tmp{j};
            end
        end
    end
end
