function DT = arbitarayPE(alpha,DT)

DT.sRPEtemp = nan(height(DT),1);
DT.sREtemp = nan(height(DT),1);

for monkey = 1:6
    for sess = 1:max(DT.session(DT.monkey==monkey))
        Q = ones(sum(DT.session==sess&DT.monkey==monkey),1)*2;
        delta = nan(sum(DT.session==sess&DT.monkey==monkey),1);
        r = DT.reward(DT.session==sess&DT.monkey==monkey);
        for t = 1:length(Q)
            delta(t) = r(t) - Q(t);
            Q(t+1) = Q(t) + alpha * delta(t);
        end
        DT.sRPEtemp(DT.session==sess&DT.monkey==monkey)=delta;
        DT.sREtemp(DT.session==sess&DT.monkey==monkey)=Q(1:end-1);
    end
end