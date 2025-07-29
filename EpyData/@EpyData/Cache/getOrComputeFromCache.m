function result = getOrComputeFromCache(obj, cacheKey, computeFcn)
% Retrieves a cached result if available, otherwise computes and stores it.
% Automatically clears the cache if size exceeds 20 GB.

maxCacheSize = 20 * 1024^3; % Size limit in bytes

% Check if result is already cached
if isfield(obj.ComputationCache, cacheKey)
    result = obj.ComputationCache.(cacheKey);
    return;
end

% Compute the result
result = computeFcn();

% Store result in cache
obj.ComputationCache.(cacheKey) = result;

% Check cache memory size
if obj.getCacheSize() > maxCacheSize
    warning('EpyData:CachePurged', 'Cache cleared (exceeded 20 GB).');
    obj.clearCache();
end
end
