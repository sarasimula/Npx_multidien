function total = getCacheSize(obj)
% Estimates the memory size of the cache (in bytes)
names = fieldnames(obj.ComputationCache);
total = 0;
for i = 1:numel(names)
    total = total + whosSize(obj.ComputationCache.(names{i}));
end
end

function s = whosSize(var)
info = whos('var');
s = info.bytes;
end