function ISstateTimestamps = IS(this, opt)
%IS Retrieve timestamps of interictal spikes (IS) for selected rats and brain states.
%
%   ISstateTimestamps = IS(obj) returns all IS timestamps for all available rats,
%   regardless of brain state.
%
%   ISstateTimestamps = IS(obj, opt) allows filtering of the IS timestamps by:
%     - rat ID(s)
%     - brain state: "sleep", "wake", or "all"
%
%   This function is useful for extracting raw IS events for downstream analyses,
%   such as temporal alignment, state-specific rate computation, or visualization.
%
%   INPUT:
%     obj     : EpyData object containing IS timestamps and brain state intervals.
%     opt     : (optional) struct with the following fields:
%         - rat   : [N x 1] double array of rat IDs. Default: all available.
%         - state : string, one of "all", "sleep", "wake". Default: "all".
%
%   OUTPUT:
%     ISstateTimestamps : Cell array of size [N x 1], containing the IS timestamps
%                         per rat, filtered by brain state.
%
%   Requires: ISTimestamps, sleepInt, wakeInt properties; Restrict function.


arguments
    this (1,1) EpyData
    opt.rat(1,:) double = this.ratID  % optionnal, by défault all the rats
    opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all" % By default, all IS regardless of the state
end

% Preallocate output cell array
ISstateTimestamps = cell(numel(opt.rat),1);
idx = 0;  % Output index counter

for j = 1:numel(opt.rat)
    IDX = idx; % Store current output index to detect matching rats

    for i = 1:numel(this.ratID)
        if this.ratID(i) == opt.rat(j)
            idx = idx + 1;
            
            % Select IS timestamps depending on the brain state
            switch opt.state
                case "all"
                    ISstateTimestamps{idx} = this.ISTimestamps{i};
                case "sleep"
                    if ~isempty(this.sleepInt{i})
                        ISstateTimestamps{idx} = Restrict(this.ISTimestamps{i},this.sleepInt{i});
                    end
                case "wake"
                    if ~isempty(this.wakeInt{i})
                        ISstateTimestamps{idx} = Restrict(this.ISTimestamps{i},this.wakeInt{i});
                    end
            end
        end
    end

     % If no match found for current rat, issue warning and leave cell empty
    if IDX == idx
        warning("No matching rat found in for rat number %s Returning empty output.",opt.rat(j));
    end

end
end

