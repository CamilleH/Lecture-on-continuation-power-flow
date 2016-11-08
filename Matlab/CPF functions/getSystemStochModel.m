function StochModel = getSystemStochModel(systemName,varargin)
StochModel = eval([systemName '_StochModel(varargin{:})']);
end