classdef MarkovModel < matlab.mixin.SetGet
    
    properties
        
        rawData;
        dataTs;
        
        % Stats data
        dataBinCounts;
        dataBinEdges;
        dataBin;
        dataBinWidth
        dataCDF;
        
        DataExpectation;
        
        randSeed;                
    end
    
    properties % Markov 1d
        markovTransMatrix = [];
        markovDimSizes = [];
        dataBinCenters = [];
        stateValidityMap = [];
    end
    

    properties
        nHistBins = 50;
        nTimesInStateIx = [];
    end
    
    methods
        
        function obj = MarkovModel()
            obj.randSeed = randi(2^31);
        end        
        
        function state = dataToState(obj,varargin)
            % Convert the state to the data
            %
            parser = inputParser();
            parser.addRequired('data',@(x)isnumeric(x)||isduration(x));
            parser.parse(varargin{:});
            p = parser.Results;
            
            state = zeros(size(obj.markovDimSizes));
            
            for idx = 1:numel(obj.markovDimSizes)
                state(idx) = interp1(obj.dataBinCenters{idx},1:numel(obj.dataBinCenters{idx}),p.data(idx),'nearest','extrap');
            end
        end 
        
        function data = stateToData(obj,varargin)
            % Convert a state into the 
            %
            parser = inputParser();
            parser.addRequired('state',@(x)isnumeric(x));
            parser.parse(varargin{:});
            p = parser.Results;
            
            data = zeros(size(obj.markovDimSizes));
            
            for idx = 1:numel(obj.markovDimSizes)
                data(idx) = interp1(1:numel(obj.dataBinCenters{idx}),obj.dataBinCenters{idx},p.state(idx),'nearest','extrap');
            end
        end
        
        function exp = getExpectation(obj,state)
            
            exp = zeros(size(obj.markovDimSizes));
            
            for idx = 1:numel(obj.markovDimSizes)
                ind = UTIL.sub2ind_vec(size(obj.DataExpectation),[state idx]);
                exp(idx) = obj.DataExpectation(ind);
            end
            
        end
        
        % Create validity map based on non-uniform time
        function computeStateValidityMap(obj,varargin)
            parser = inputParser();
            parser.addRequired('maxStepSize',@(x)isnumeric(x)||isduration(x));
            parser.addParameter('finiteCheck',true,@islogical);
            parser.parse(varargin{:});
            p = parser.Results;
            
            
            if(istimetable(obj.rawData))
                time = obj.rawData.Properties.RowTimes;
                data = table2array(obj.rawData);
                assert(isduration(p.maxStepSize),'MarkovModel::computeStateValidityMatrix - Max Step Size must be a duration type for timetable data');
                minBetweenDataSamp = diff(time);
                
            elseif(isa(obj.rawData,'timeseries'))
                time = obj.rawData.Time;
                data = obj.rawData.Data;
                if(isduration(p.maxStepSize))
                    warning('MarkovModel::computeStateValidityMatrix - Max step size provided as duration but time series data has no datetime type. Assuming datenum');
                    minBetweenDataSamp = diff(datetime(time,'ConvertFrom','datenum'));
                else
                    minBetweenDataSamp = diff(time);
                end
                
            elseif(isnumeric(obj.rawData))
                error('MarkovModel::computerStateValidityMap - Cannot genrate validity map with numeric input');
            end
            
            % Since the data is not necessarily regularly sampled and there
            % are large gaps, we need to exclude those
            obj.stateValidityMap = [minBetweenDataSamp <= p.maxStepSize; false];
            
            if(p.finiteCheck)
                obj.stateValidityMap = obj.stateValidityMap & all(isfinite(data),2);
            end
        
            
        end
        
        % Process inputs to compute transition matrix
        function computeNdMarkovMatrix(obj,varargin)
            parser = inputParser();
            parser.addParameter('forceRegen',false,@islogical);
            parser.parse(varargin{:});
            p = parser.Results;
            
            if isempty(obj.dataBin)
                error('MarkovModel::compute2dMarkovMatrix - Bins not initialized, run createStateBins()');
            end
            
            % Generate the markov transition matrix for the wind generation
            if isempty(obj.markovTransMatrix) || p.forceRegen
                [obj.markovTransMatrix, obj.nTimesInStateIx] = obj.genNdMarkovMatrix(obj.dataBin,obj.stateValidityMap,obj.markovDimSizes);
            end
            
        end
        
        % Calculate the expectation matrix
        function generateExpectationMatrix(obj,varargin)
            parser = inputParser();
            %parser.addParameter('forceRegen',false,@islogical);
            parser.parse(varargin{:});
            p = parser.Results;
            
            % Pre-calculate the expectation vector based on the markov
            % transition matrix by sum(P_ij * state_ij)
            obj.DataExpectation = zeros([obj.markovDimSizes, numel(obj.markovDimSizes)]);
            
            % Get the elements of the transition matrix - Since its sparse
            % we need to use find to return the i,j index and the value (s)
            [i,j,s] = find(obj.markovTransMatrix);
            
            % Loop over unique starting states
            ui = unique(i);
            for idx = 1:numel(ui)
                
                i_idx = find(i == ui(idx));
                
                % Extract the transition probability
                probs = s(i_idx);
                
                destStates = j(i_idx);
                
                % Decode the linear state into vector components of each dimension
                decodeState = UTIL.ind2sub_vec(obj.markovDimSizes,destStates);
                
                % Populate the expectation matrix by multiplying each
                % terminal node value by the probability and summing
                for jdx = 1:numel(obj.markovDimSizes)
                    tmp = obj.dataBinCenters{jdx}(decodeState(:,jdx));
                    obj.DataExpectation(ui(idx) + (jdx-1) * prod(obj.markovDimSizes)) = sum(tmp .* probs);
                end
            end
        end
        
        % Calculate the expectation matrix
        function generateFullExpectationMatrix(obj,varargin)
            parser = inputParser();
            %parser.addParameter('forceRegen',false,@islogical);
            parser.parse(varargin{:});
            p = parser.Results;
            
            % Pre-calculate the expectation vector based on the markov
            % transition matrix by sum(P_ij * state_ij)
            obj.DataExpectation = zeros([obj.markovDimSizes, numel(obj.markovDimSizes)]);
            
            % Get the elements of the transition matrix - Since its sparse
            % we need to use find to return the i,j index and the value (s)
            [i,j,s] = find(obj.markovTransMatrix);
            
            % Loop over unique starting states
            ui = unique(i);
            
            % Create a matrix of NaNs with the index of the populated
            % states found
            ui_full = nan(1,size(obj.markovTransMatrix,1));
            ui_full(i) = i;
            
            % If it is 1-dimension reshape won't work, so test it
            if numel(obj.markovDimSizes) > 1
                ui_full_nd = reshape(ui_full,obj.markovDimSizes);
            else
                ui_full_nd = ui_full;
            end
            
            % Replace NaNs with nearest neighbor
            nonnan = ~isnan(ui_full_nd);
            vvalues = ui_full_nd(nonnan);
            [~,y] = min(abs(bsxfun(@minus,1:numel(ui_full_nd),find(nonnan))));
            %ui_full_nd(~nonnan) = vvalues(y(~nonnan));
            
            for idx = 1:size(obj.markovTransMatrix,1)
                
                
                if any(ui == idx) % If this state was found somewhere
                    i_idx = find(i == idx);
                    
                else % If not, find the nearest neighbor
                    % This may not be a dimensionally valid assumption
                    i_idx = ui_full_nd(idx);
                end
                
                if isnan(i_idx); i_idx = 1; end
                
                
                % Extract the transition probability
                if (~isempty(s))
                    probs = s(i_idx);
                else
                    probs = 0;
                end

                if (~isempty(j))
                    destStates = j(i_idx);
                else
                    destStates = 0;
                end

                % Decode the linear state into vector components of each dimension
                decodeState = UTIL.ind2sub_vec(obj.markovDimSizes,destStates);
                
                % Populate the expectation matrix by multiplying each
                % terminal node value by the probability and summing
                for jdx = 1:numel(obj.markovDimSizes)
                    tmp = obj.dataBinCenters{jdx}(decodeState(:,jdx));
                    obj.DataExpectation(idx + (jdx-1) * prod(obj.markovDimSizes)) = sum(tmp .* probs);
                end
            end
        end
                        
        % Treat raw data as pre-binned
        function useDataAsBins(obj,varargin)
            parser = inputParser();
            parser.addOptional('dataBins',[],@isnumeric);
            parser.parse(varargin{:});
            p = parser.Results;
            
            if(istimetable(obj.rawData))
                data = table2array(obj.rawData);
            elseif(isa(x,'timeseries'))
                data = obj.rawData.Data;
            elseif(isnumeric(obj.rawData))
                data = obj.rawData;
            end
            
            obj.dataBinCounts = {};
            obj.dataBinEdges = {};
            obj.dataCDF = {};
            obj.dataBin = [];
            obj.dataBinCenters = [];
            for sdx = 1:size(data,2)
                
                tmp = data(:,sdx);
                
                obj.dataBinCounts{sdx} = unique(data(:,sdx));
                obj.dataBinEdges = unique(data(:,sdx));
                obj.dataCDF = {};
                
            end
            
            error('MarkovModel::useDataAsBins - Not implemented yet');
            
        end
        
        % Bin the data into discrete states
        function createStateBins(obj,varargin)
            parser = inputParser();
            parser.addOptional('histBins',[],@isnumeric);
            parser.parse(varargin{:});
            p = parser.Results;
            
            if(istimetable(obj.rawData))
                data = table2array(obj.rawData);
            elseif(isa(x,'timeseries'))
                data = obj.rawData.Data;
            elseif(isnumeric(obj.rawData))
                data = obj.rawData;
            end
            
            % Normalize bin size
            if(size(data,2) > 1)
                if(numel(p.histBins) == 1)
                    p.histBins = repmat(p.histBins,1,size(data,2));
                elseif(numel(p.histBins) > size(data,2))
                    error('MarkovModel::createStateBins - Number of hist bins provided may not be greater than number of data columns');
                end
                
            end
            
            obj.dataBinCounts = {};
            obj.dataBinEdges = {};
            obj.dataCDF = {};
            obj.dataBin = [];
            obj.dataBinCenters = [];
            for sdx = 1:size(data,2)
                
                % Calculate statistics
                if(isempty(p.histBins) || p.histBins(sdx) == -1)
                    [tmpBinCounts, tmpEdges, tmpBins] = histcounts(data(:,sdx), 'Normalization', 'probability');
                else
                    [tmpBinCounts, tmpEdges, tmpBins] = histcounts(data(:,sdx), p.histBins(sdx), 'Normalization', 'probability');
                end
                
                if(isrow(tmpBinCounts))
                    tmpBinCounts = reshape(tmpBinCounts,[],1);
                end
                obj.dataBinCounts{end+1} = tmpBinCounts;
            
                if(isrow(tmpEdges))
                    tmpEdges = reshape(tmpEdges,[],1);
                end
                obj.dataBinEdges{end+1} = tmpEdges;
                
                % Calculate center of each data bin
                obj.dataBinCenters{end+1} = tmpEdges(1:end-1) + diff(tmpEdges([1,2]))/2;
                
                
                datCdf = cumtrapz(tmpBinCounts);
                % Make sure CDF values are unique
                if numel(unique(datCdf)) ~= numel(datCdf)
                    % TODO: Fix this to be a more scalable approach
                    for idx = 1:numel(datCdf)
                        datCdf(idx) = datCdf(idx)+idx*0.00001;
                    end
                end
                obj.dataCDF{end+1} = datCdf;
                
                if(isrow(tmpBins))
                    tmpBins = reshape(tmpBins,[],1);
                end
                
                % Remove any invalid entries like NaNs
                tmpBins = tmpBins(not(tmpBins == 0));
                obj.dataBin = [obj.dataBin tmpBins];
                
            end           
            
            % Set the dimensions of the state matrix
            obj.markovDimSizes = cellfun(@numel,obj.dataBinCounts);
            
        end
        
        function removeTerminalStates(obj,varargin)
            parser = inputParser();
            %parser.addOptional('histBins',[],@isnumeric);
            parser.parse(varargin{:});
            p = parser.Results;
            
            iter = 0;
            while(true)
                
                %obj.plot2dStates('fignum',iter+1);
                
                iter = iter + 1;
                % Get the elements of the transition matrix - Since its sparse
                % we need to use find to return the i,j index and the value (s)
                [i,j,~] = find(obj.markovTransMatrix);
                
                % Find source and destination elements that dont have
                % counterparts
                orphanI = setdiff(unique(i),unique(j));
                
                % Dead ends - Destination states in J that are not present in I
                orphanJ = setdiff(unique(j),unique(i));
                
                if numel(orphanJ) == 0 && numel(orphanI) == 0
                    break;
                end
                
                fprintf(1,'Iteration %d, %d terminal states found\n',iter,numel(orphanJ));
                
                for odx = 1:numel(orphanJ)
                    tmp = UTIL.ind2sub_vec(obj.markovDimSizes, orphanJ(odx));
                    
                    % Look backwards
                    obj.findLastNodeWithMultiplePaths(orphanJ(odx),0);
                    % Look forwards
                    obj.findLastNodeWithMultiplePaths(orphanJ(odx),1);
                    
                end
                
                fprintf(1,'Iteration %d, %d origin states found\n',iter,numel(orphanI));
                for odx = 1:numel(orphanI)
                    tmp = UTIL.ind2sub_vec(obj.markovDimSizes, orphanI(odx));
                    
                    % Look forwards
                    obj.findLastNodeWithMultiplePaths(orphanI(odx),1);
                    % Look backwards
                    obj.findLastNodeWithMultiplePaths(orphanI(odx),0);
                end
                
            end

        end
        
        function findLastNodeWithMultiplePaths(obj,linIx,dir)
            
            orgNode = UTIL.ind2sub_vec(obj.markovDimSizes, linIx);
            
            if dir == 0 % backwards
                
                % Find nodes this goes to (should be none)
                [~,j,~] = find(obj.markovTransMatrix(linIx,:));
                if numel(j) > 1
                    warning('MarkovModel::findLastNodeWithMultiplePaths - Attempting to remove node with valid paths, why?? Returned');
                   return; % We should never find more than one outgoing path here. If so, return
                end
                
                % Find all nodes where this is the next node
                [i,~,~] = find(obj.markovTransMatrix(:,linIx));
                
                prevNodesFound = numel(i);
                
                if prevNodesFound == 0 % We were on a line, remove it all
                    obj.markovTransMatrix(linIx,:) = 0;
                    fprintf(1,'Removing terminal/origin node %s\n',mat2str(orgNode));
                else
                
                    % If there are multiple org nodes, process both
                    for nix = 1:numel(i)
                        % Check to see if nodes in i have other outgoing paths
                        %  If not, remove node and rescale probabilities at i
                        [~,j,~] = find(obj.markovTransMatrix(i(nix),:));

                        altPathsFound = numel(j) - 1;

                        % If no alternate paths, keep searching backwards
                        if altPathsFound == 0

                            fprintf(1,'Removing terminal node %s and continuing\n',mat2str(orgNode));
                            obj.findLastNodeWithMultiplePaths(i(nix),dir)
                            % If we got here, it's also a terminal node, so remove it
                            obj.markovTransMatrix(linIx,:) = 0;
                        elseif altPathsFound < 0 
                            % shouldnt hit this but just in case
                            error('MarkovModel::findLastNodeWithMultiplePaths - Alternate paths found < 0. Should be invalid');

                        else % If there are other outgoing paths, remove node I and rescale probs
                            fprintf(1,'Removing terminal node %s and rescaling incoming probabilities\n',mat2str(orgNode));
                            obj.removeTerminalNodeAndRescaleProbabilities(i(nix),linIx);
                        end
                    end
                
                end
                
                
            else % forwards
                [~,j,s] = find(obj.markovTransMatrix(linIx,:));
                
                nextNodesFound = numel(j);
                if nextNodesFound == 1
                    obj.findLastNodeWithMultiplePaths(j,dir)
                elseif nextNodesFound > 1
                    % We found a node with multiple other next nodes. Stop
                else
                    % We reached a dead end
                end
                
                if ~isempty(s) % If this node is populated
                    % No matter how we get here, remove the node
                    fprintf(1,'Removing terminal node %s and continuing\n',mat2str(orgNode));
                    obj.markovTransMatrix(linIx,:) = 0;
                end
            end
            
        end
        
        function removeTerminalNodeAndRescaleProbabilities(obj,srcNode,nodeToRemove)
            [~,j,s] = find(obj.markovTransMatrix(srcNode,:));
            
            % Find previous elements to be removed
            ix = (j == nodeToRemove);
            
            if numel(ix) == 1 % This might happen ?
                obj.markovTransMatrix(srcNode,nodeToRemove) = 0;
            else
                % Rescale the probability and evenly apply it to remaining elements
                s(~ix) = s(~ix) + (sum(s(ix)) / sum(~ix));

                % Update transition matrix and pop off terminal node
                obj.markovTransMatrix(srcNode,j(~ix)) = s(~ix);
                obj.markovTransMatrix(srcNode,j(ix)) = 0;
            end
            
        end
        
    end
    
    methods % Data generation
       
        function dataSamp = genSampleData(obj,varargin)
            parser = inputParser();
            parser.addRequired('numSamples',@(x)isnumeric(x)&&x>0);
            parser.addOptional('StartingTime',datetime('1/1/2020 12:00:00 AM','InputFormat','MM/dd/uuuu hh:mm:ss aa'));
            parser.addOptional('StepSize_min',10);
            parser.addOptional('initSample',[]);
            parser.addParameter('initialMarkovStates',[],@(x)isnumeric(x) && all(x > 0));
            parser.parse(varargin{:});
            p = parser.Results;
            
            % Set a seed so the results are repeatable
            if(~isempty(obj.randSeed))
                rng(obj.randSeed);
            end
                        
            % If not given, assume starting state is the mean
            if isempty(p.initialMarkovStates)
                linStates = UTIL.sub2ind_vec(obj.markovDimSizes,obj.dataBin);
                iStates = mode(linStates);
                iStates = UTIL.ind2sub_vec(obj.markovDimSizes,iStates);
            else
                assert(numel(p.initialMarkovStates) == numel(obj.markovDimSizes),'MarkovModel::genSampleData - Initial states must be provided for each dimension')
                iStates = p.initialMarkovStates;
            end
            
            % Map the initial state to the linear index
            linState = UTIL.sub2ind_vec(obj.markovDimSizes, iStates);
            
            % Preallocate output matrix
            nStates = numel(obj.markovDimSizes);
            dataSamp = zeros(p.numSamples,nStates);
            
            for sdx = 1:p.numSamples
                
                % Decode the linear state into i,j,... states
                tmp = UTIL.ind2sub_vec(obj.markovDimSizes, linState);
                
                for ddx = 1:nStates
                   dataSamp(sdx,ddx) = obj.dataBinCenters{ddx}(tmp(ddx));
                end
                
                % Extract transition probabilities from the transition and
                % compute the CDF to inverse sample against
                [~,nextIx,probs] = find(obj.markovTransMatrix(linState,:));
                probs = cumsum(probs);
                
                % Find next sample
                nextSample = find(rand() < probs,1,'first');
                
                if(~isempty(nextSample))
                    linState = nextIx(nextSample);
                else
                    warning('MarkovModel::genSampleData - Entered terminal state. Iteration %d, source state %s',sdx,mat2str(tmp))
                    break;
                end
            end
            
            % Reset rng if needed
            if(~isempty(obj.randSeed))
                rng('default');
            end
        end 
        
    end
    
    methods % Data ingestion
        
        
        function importTimeseries(obj,varargin)
            parser = inputParser();
            parser.addRequired('dataTs',@(x)isa(x,'timeseries'));
            parser.addOptional('varName','',@ischar);
            parser.parse(varargin{:});
            p = parser.Results;
            
            warning('MarkovModel::importTimeseries - May not be fully implemented');
            
            obj.rawData = p.dataTs;
            
            % Default state validity to true
            obj.stateValidityMap = true(size(obj.rawData,1),1);
        end
        
        function importTimeTable(obj,varargin)
            parser = inputParser();
            parser.addRequired('dataTs',@istimetable);
            parser.addRequired('dataVarName',@(x)ischar(x)||iscellstr(x)||isstring(x));
            parser.addParameter('sampleTime',[],@isduration);
            parser.parse(varargin{:});
            p = parser.Results;
            
            if ischar(p.dataVarName)
                p.dataVarName = {p.dataVarName};
            end
            
            assert(any(cellfun(@(x)any(strcmpi(x,p.dataTs.Properties.VariableNames)),p.dataVarName)),'MarkovModel::importTimeTable - Data column not found in time table');
            
            % Find columns and remove elements not found
            
            idxs = false(size(p.dataTs.Properties.VariableNames));
            for idx = p.dataVarName
                idxs = idxs | strcmpi(idx,p.dataTs.Properties.VariableNames);
            end
            
            obj.rawData = p.dataTs(:,idxs);
            
            % Default state validity to true
            obj.stateValidityMap = true(size(obj.rawData,1),1);
            
        end
        
        function importArray(obj,varargin)
            parser = inputParser();
            parser.addRequired('dataTs');
            parser.addOptional('varName','',@ischar);
            parser.addParameter('scaleWind',1,@isnumeric);
            parser.parse(varargin{:});
            p = parser.Results;
            
            assert(false,'MarkovModel::importArray - Not implemented');
        end
        
    end
    
    properties
        
        activeLines = matlab.graphics.primitive.Line.empty;
        
    end
    
    methods % Plotting methods
                
        function plotExpectationMatrix(obj,varargin)
            parser = inputParser();
            parser.addOptional('fignum',380,@isnumeric);
            parser.parse(varargin{:});
            p = parser.Results;
            
            figure(p.fignum)
            ax_im = subplot(1,1,1);
            if( numel(obj.markovDimSizes) == 3)
                error('MarkovModel::plotExpectationMatrix - Cannot plot expectation matrix of > 2 dimensions');
            elseif( numel(obj.markovDimSizes) == 2)
                surf(ax_im,obj.DataExpectation(:,:,1))
                colorbar
            elseif( numel(obj.markovDimSizes) == 1)
                plot(ax_im,obj.DataExpectation,'bo','LineStyle','-')
            end
            
            title(ax_im,'Markov Expectation Matrix');
            xlabel(ax_im,'State i');
            ylabel(ax_im,'State j');
            
        end

        %% TODO: Implement this
        % function plotOneHistogram(obj,varargin)
        %     parser = inputParser();
        %     parser.addOptional('fignum',380,@isnumeric);
        %     parser.parse(varargin{:});
        %     p = parser.Results;
        %     
        %     figure(p.fignum)
        %     ax_im = subplot(1,1,1);
        %     if( numel(obj.markovDimSizes) == 3)
        %         error('MarkovModel::plotExpectationMatrix - Cannot plot expectation matrix of > 2 dimensions');
        %     elseif( numel(obj.markovDimSizes) == 2)
        %         surf(ax_im,obj.DataExpectation(:,:,1))
        %         colorbar
        %     elseif( numel(obj.markovDimSizes) == 1)
        %         plot(ax_im,obj.DataExpectation,'bo','LineStyle','-')
        %     end
        %     
        %     title(ax_im,'Markov Expectation Matrix');
        %     xlabel(ax_im,'State i');
        %     ylabel(ax_im,'State j');
        %     
        % end
        
        function addLinesToStates(obj,src,evtData)
            
            if ~isempty(obj.activeLines)
                for l = obj.activeLines
                    delete(obj.activeLines);
                    obj.activeLines = matlab.graphics.primitive.Line.empty;
                end
            end
            
            ndim = numel(obj.markovDimSizes);
            
            if(ndim == 3)
                pointClicked = fliplr(evtData.IntersectionPoint);
            elseif(ndim ==2)
                pointClicked = fliplr(evtData.IntersectionPoint(1:2));
            elseif(ndim == 1)
                pointClicked = evtData.IntersectionPoint(1);
%                 error('MarkovModel::addLinesToStates - Cant do 1 dimension yet');
            else
                error('MarkovModel::addLinesToStates - Invalid dimensions');
            end
            
            srcState = UTIL.sub2ind_vec(obj.markovDimSizes,pointClicked);
            
            [~,j,s] = find(obj.markovTransMatrix(srcState,:));
            
            % Plot outgoing lines
            for ldx = 1:numel(j)
                start_end = UTIL.ind2sub_vec(obj.markovDimSizes, [srcState;j(ldx)]);
                
                if( size(start_end,2) == 3)
                    obj.activeLines(end+1) = line(src.Parent,start_end(:,3),start_end(:,2),start_end(:,1),'LineWidth',max(numel(s)*s(ldx),1),'Color','g');
                elseif( size(start_end,2) == 2)
                    obj.activeLines(end+1) = line(src.Parent,start_end(:,2),start_end(:,1),'LineWidth',max(numel(s)*s(ldx),1),'Color','g');
                elseif( size(start_end,2) == 1)
                    obj.activeLines(end+1) = line(src.Parent,start_end(:,1),start_end(:,1),'LineWidth',max(numel(s)*s(ldx),1),'Color','g');
                end
            end
            
            [i,~,s] = find(obj.markovTransMatrix(:,srcState));
            
            % Plot incoming lines
            for ldx = 1:numel(i)
                start_end = UTIL.ind2sub_vec(obj.markovDimSizes, [srcState;i(ldx)]);
                
                if( size(start_end,2) == 3)
                    obj.activeLines(end+1) = line(src.Parent,start_end(:,3),start_end(:,2),start_end(:,1),'LineWidth',max(numel(s)*s(ldx),1),'Color','r');
                elseif( size(start_end,2) == 2)
                    obj.activeLines(end+1) = line(src.Parent,start_end(:,2),start_end(:,1),'LineWidth',max(numel(s)*s(ldx),1),'Color','r');
                elseif( size(start_end,2) == 1)
                    obj.activeLines(end+1) = line(src.Parent,start_end(:,1),start_end(:,1),'LineWidth',max(numel(s)*s(ldx),1),'Color','r');
                end
            end
            
        end
        
        function plot2dStates(obj,varargin)
            
                parser = inputParser();
                parser.addOptional('fignum',399,@isnumeric);
                parser.parse(varargin{:});
                p = parser.Results;
            
                fig_h = figure(p.fignum);
                clf(fig_h);
                
                ax_h = subplot(1,1,1);
                hold(ax_h,'on');
                
                % Find all transitions in the matrix
                [i,j,s] = find(obj.markovTransMatrix);
                
                i_uniq = unique(i);
                j_uniq = unique(j);
                subIxs = repmat(zeros(length(i_uniq),1),1,numel(obj.markovDimSizes)+1);
                subJxs = repmat(zeros(length(j_uniq),1),1,numel(obj.markovDimSizes)+1);
                
                for idx = 1:length(i_uniq)
                    subIxs(idx,1:end-1) = UTIL.ind2sub_vec(obj.markovDimSizes, i_uniq(idx));
                    subIxs(idx,end) = numel(find(obj.markovTransMatrix(i_uniq(idx),:)));
                end
                
                for jdx = 1:length(j_uniq)
                    subJxs(jdx,1:end-1) = UTIL.ind2sub_vec(obj.markovDimSizes, j_uniq(jdx));
                    subJxs(jdx,end) = numel(find(obj.markovTransMatrix(j_uniq(jdx),:)));
                end
                
                % Find nodes with no termination or origination
                [i_notinJ,i_notinJ_ix] = setdiff(i_uniq,j_uniq,'stable');
                [j_notinI,j_notinI_ix] = setdiff(j_uniq,i_uniq,'stable');
                
                subInJxs = repmat(zeros(length(i_notinJ),1),1,numel(obj.markovDimSizes)+1);
                subJnIxs = repmat(zeros(length(j_notinI),1),1,numel(obj.markovDimSizes)+1);
                
                for idx = 1:length(i_notinJ)
                    subInJxs(idx,1:end-1) = UTIL.ind2sub_vec(obj.markovDimSizes, i_notinJ(idx));
                    subInJxs(idx,end) = numel(find(obj.markovTransMatrix(i_notinJ(idx),:)));
                end
                subInJxs(:,end) = 50;
                
                for jdx = 1:length(j_notinI)
                    subJnIxs(jdx,1:end-1) = UTIL.ind2sub_vec(obj.markovDimSizes, j_notinI(jdx));
                    subJnIxs(jdx,end) = numel(find(obj.markovTransMatrix(j_notinI(jdx),:)));
                end
                subJnIxs(:,end) = 60;
                
                %
                maxMarkerSize = 50;
                minMarkerSize = 10;
                tmpI = obj.normalizeMarkerSize(subIxs(:,end),minMarkerSize,maxMarkerSize);
                tmpJ = obj.normalizeMarkerSize(subJxs(:,end),minMarkerSize,maxMarkerSize);
                
                if numel(obj.markovDimSizes) == 2
                
                    sNodes = scatter(ax_h,subIxs(1:end,2),subIxs(1:end,1),floor(tmpI)+3,'bo');
                    dNodes = scatter(ax_h,subJxs(1:end,2),subJxs(1:end,1),floor(tmpJ),'rs','MarkerFaceColor','r');
                    
                    if ~isempty(subJnIxs)
                        jNiNodes = scatter(ax_h,subInJxs(1:end,2),subInJxs(1:end,1),subInJxs(:,end),'yo','MarkerFaceColor','y','MarkerEdgeColor','k');
                    end
                    
                    if ~isempty(subJnIxs)
                        iNjNodes = scatter(ax_h,subJnIxs(1:end,2),subJnIxs(1:end,1),subJnIxs(:,end),'md','MarkerFaceColor','m','MarkerEdgeColor','k');
                    end
                    
                elseif numel(obj.markovDimSizes) == 3
                    sNodes = scatter3(ax_h,subIxs(1:end,3),subIxs(1:end,2),subIxs(1:end,1),floor(tmpI)+3,'bo');
                    dNodes = scatter3(ax_h,subJxs(1:end,3),subJxs(1:end,2),subJxs(1:end,1),floor(tmpJ),'rs','MarkerFaceColor','r');

                    if ~isempty(subInJxs)
                        jNiNodes = scatter3(ax_h,subInJxs(1:end,3),subInJxs(1:end,2),subInJxs(1:end,1),subInJxs(:,end),'yo','MarkerFaceColor','y','MarkerEdgeColor','k');
                    end
                    
                    if ~isempty(subJnIxs)
                        iNjNodes = scatter3(ax_h,subJnIxs(1:end,3),subJnIxs(1:end,2),subJnIxs(1:end,1),subJnIxs(:,end),'md','MarkerFaceColor','m','MarkerEdgeColor','k');
                    end
                    
                elseif numel(obj.markovDimSizes) == 1
                    sNodes = scatter(ax_h,subIxs(1:end,1),subIxs(1:end,1),floor(tmpI)+3,'bo');
                    dNodes = scatter(ax_h,subJxs(1:end,1),subJxs(1:end,1),floor(tmpJ),'rs','MarkerFaceColor','r');

                    if ~isempty(subInJxs)
                        jNiNodes = scatter(ax_h,subInJxs(1:end,1),subInJxs(1:end,1),subInJxs(:,end),'yo','MarkerFaceColor','y','MarkerEdgeColor','k');
                    end
                    
                    if ~isempty(subJnIxs)
                        iNjNodes = scatter(ax_h,subJnIxs(1:end,1),subJnIxs(1:end,1),subJnIxs(:,end),'md','MarkerFaceColor','m','MarkerEdgeColor','k');
                    end
                else
                    error('MarkovModel::PlotNdStates - Cannot plot in %d dimensions. Please download additional spacetime',numel(obj.markovDimSizes));
                end
                
                sNodes.ButtonDownFcn = @obj.addLinesToStates;
                dNodes.ButtonDownFcn = @obj.addLinesToStates;
                jNiNodes.ButtonDownFcn = @obj.addLinesToStates;
                iNjNodes.ButtonDownFcn = @obj.addLinesToStates;
                
                title(ax_h,'Markov State Matrix');
                xlabel(ax_h,'State i');
                ylabel(ax_h,'State j');
        end
        
        function siz = normalizeMarkerSize(obj,vec,minMarkerSize,maxMarkerSize)
            siz = vec - min(vec);
            siz = siz/max(siz);
            siz = (siz * (maxMarkerSize-minMarkerSize)) + minMarkerSize;
        end
        
        function plotStateActivity(obj)
            % Plot the heatmap of state presence
            [x,y] = meshgrid(1:obj.markovDimSizes(2),1:obj.markovDimSizes(1));
            figure(22)
            dataVals = reshape(obj.nTimesInStateIx,obj.markovDimSizes(1),obj.markovDimSizes(2));
            loggedVals = log(1.1-dataVals/max(max(dataVals)));
            loggedVals(loggedVals>-0.2) = NaN;
            %imagesc((1:obj.markovDimSizes(2))/60,obj.dataBinEdges(1:obj.markovDimSizes(1)),loggedVals)
            surf(x,obj.dataBinEdges(y),dataVals);
            
            colorbar
        end
        
        function plotBinScatter(obj)

                figure(24)
                clf
                h = binscatter(gid(timeState)/60+1,obj.dataBinEdges(obj.dataBin)');
                colormap('parula');
                ax = gca;
                %ax.ColorScale = 'log';
                
                % Calculate the center of mass of each column
                Ycenters = movmean(h.YBinEdges,2);
                Ycenters = Ycenters(2:end);
                
                Xcenters = movmean(h.XBinEdges,2);
                Xcenters = Xcenters(2:end);
                
                avgLineY = [];
                for hdx = 1:size(h.Values,1)
                    nonNanVals = h.Values(hdx,:)
                    fidx = isfinite(nonNanVals);
                    nonNanVals = nonNanVals(fidx);
                    avgLineY(hdx) = sum(nonNanVals.*Ycenters(fidx)) / sum(nonNanVals);
                end
                hold on;
                plot(Xcenters,avgLineY,'r-','LineWidth',4,'DisplayName','Mean State Occurrences');
                
                title('State Distribution and Mean');
        end
        
        function plotTransitionMatrix(obj,varargin)
            
                parser = inputParser();
                parser.addOptional('transitionMatrix',obj.markovTransMatrix);
                parser.parse(varargin{:});
                p = parser.Results;
            
                figure(453);
                clf;
                ax_im = subplot(1,1,1);
                imagesc(ax_im,p.transitionMatrix);
                colorbar(ax_im);
                colormap('jet');
                title(ax_im,'Markov Transition Matrix');
                xlabel(ax_im,'State i');
                ylabel(ax_im,'State j');
            
        end
        
        function plotTransitionMatrixDensity(obj)
            figure(23)
            spy(obj.markovTransMatrix);
            title('Markov Transition Matrix Density');
            xlabel('State i');
            ylabel('State j');
        end
        
    end
    
    
    methods (Static)
        
        function [nodeTxMatrix,nTimesInStateIx] = genNdMarkovMatrix(stateVector,stateValidityMap,state_bin_size)
            
            assert(size(stateVector,1) > size(stateVector,2),'genNdMarkovMatrix - State vector must be column vectors');
            
            if isempty(stateValidityMap)
                stateValidityMap = true(length(stateVector));
            end
            assert(all(length(stateValidityMap) == length(stateVector)),'genNdMarkovMatrix - Validity map and state vector must be equal length');
            
            % For each dimension, determine the size of the state
            dimSize = state_bin_size;
%             dimSize = zeros(1,size(stateVector,2));
            uniqueStateVals = {};
            for idx = 1:numel(dimSize)
                
                % Pack the unique values for each dimension so we don't have to
                % keep searching for it
%                 utmp = unique(stateVector(:,idx));
                utmp = 1:dimSize(idx);
                uniqueStateVals{idx} = utmp;
                
                % Store the size in a separate array
%                 dimSize(idx) = size(utmp,1);
            end
            
            % Calculate number of unique state combinations
            nStateComb = prod(dimSize);
            
            % Collect count of times in each state comb
            nTimesInStateIx = zeros(1,nStateComb);
            
            % Find allocations
            for sdx = 1:nStateComb
                % Find each instance in time where the system was in state (i,j)
                dimVec = UTIL.ind2sub_vec(dimSize,sdx);
                
                % Find occurrences of that state combination
                iIdxs = all((stateVector == dimVec),2);
                nTimesInStateIx(sdx) = sum(iIdxs);
            end
                        
            % Look at the number of samples of this bin relative to others.
            % If there are significantly more, provide an allocation and
            % randomly discard any beyond that. This may or may not be used
            % depending on the density of each state occurrence
            allocation = floor(std(nTimesInStateIx)*3);
            
            % Allocate transition matrix
            % Give it an additional dimension to store the TO state
            %   -- It may need an extra [n1,n2,nN,... [1,nN-1]] dimensions?
            nodeTxMatrix = sparse(nStateComb,nStateComb);
            
            % nodeTxMatrix is indexed as nodeTxMatrix( S1_k, S2_k, pS1_k+1,
            % pS2_k+1
            for sdx = 1:nStateComb
                
                % If that state combination is found
                if nTimesInStateIx(sdx) > 0
                    
                    % Find the row,col index based on the linear index and dimension
                    sIdx = UTIL.ind2sub_vec(dimSize,sdx);
                    
                    fprintf(1,'State %s, %d times in state.\n',mat2str(sIdx),nTimesInStateIx(sdx));
                    
                    stateVals = zeros(size(dimVec));
                    for ddx = 1:numel(dimSize)
                        stateVals(ddx) = uniqueStateVals{ddx}(sIdx(ddx));
                    end
                    
                    % Find occurrences of that state combination
                    iIdxs = all((stateVector == stateVals),2);
                    assert(sum(iIdxs) > 0,'We thought something was in this state but nothing was found.');
                    
                    % Filter out invalid states
                    iIdxs = iIdxs & stateValidityMap;
                    
                    % If there are no valid states, skip this iteration
                    if(sum(iIdxs) == 0)
                        continue;
                    end
                    
                    % If the state is overallocated, reduce it to a 
                    %if(nTimesInStateIx(sdx) > allocation)
                    %    fprintf(1,' * State %s overallocated - reducing samples to %d.\n',mat2str(sIdx),allocation);
                    %    indxNum = randsample(find(iIdxs),nTimesInStateIx(sdx)-allocation);
                    %    iIdxs(indxNum) = false;
                    %end
                    
                    % Shift indicies by 1 to find the k+1 state
                    iPlusIdxs = [false; iIdxs(1:end-1)];
                    
                    % Assume the next state occurrences are << than the
                    % total number of samples in our input dataset. Find
                    % the unique values
                    nextStates = unique(stateVector(iPlusIdxs,:),'rows');
                    nextCounts = groupcounts(stateVector(iPlusIdxs,:));
                    nextProbs = nextCounts / sum(nextCounts);
                    
                    % Add each element to the transition matrix
                    for ndx = 1:size(nextStates,1)

                        destInd = nextStates;
                        destSub = UTIL.sub2ind_vec(dimSize,nextStates(ndx,:));
                        
                        nodeTxMatrix(sdx,destSub) = nextProbs(ndx);
                        
                    end
                end
            end
        end
        
        function P = genMarkovMatrix(binArr,nBins,validityMap)
            P = zeros(nBins,nBins);
            
            
            % Find allocations
            for idx = 1:nBins
                % Find each instance in time where the system was in state i
                iIdxs = (binArr == idx);
                nTimesInStateI(idx) = sum(iIdxs);
            end
            
            % Look at the number of samples of this bin relative to others.
            % If there are significantly more, provide an allocation and
            % randomly discard any beyond that.
            allocation = floor(std(nTimesInStateI)*3);
            
            for idx = 1:nBins
                % Find each instance in time where the system was in state i
                iIdxs = (binArr == idx);
                nTimesInStateI = sum(iIdxs);
                
                % Remove instances where the next sample is invalid
                iIdxs = iIdxs & [validityMap(2:end); false];                
                nValidTimesInState_I = sum(iIdxs);
                
                if(nValidTimesInState_I > allocation)
                    indxNum = randsample(find(iIdxs),nValidTimesInState_I-allocation);
                    iIdxs(indxNum) = false;
                end
                
                iPlusIdxs = [false; iIdxs(1:end-1)];
                
                fprintf(1,'State %d, %d times in state, %0.2f%% valid.\n',idx,nTimesInStateI,nValidTimesInState_I/nTimesInStateI*100);
                
                for jdx = 1:nBins
                    
                    % Look at bin i and determine the probability it will transition to bin j to get
                    % P_ij
                    nTimesItransToJ = sum(binArr(iPlusIdxs) == jdx);
                    P(idx,jdx) = nTimesItransToJ/nValidTimesInState_I;
                end
            end
        end
        
        function [fig_h, ax_ws] = plotDataWithDistribution(varargin)
            parser = inputParser();
            parser.addRequired('dataSamp',@isnumeric);
            parser.addOptional('time',[],@(x)isnumeric(x)||isdatetime(x)||isduration(x));
            parser.addParameter('title','Data',@ischar);
            parser.addParameter('plotName','',@ischar);
            parser.addParameter('dataName','Value',@ischar);
            parser.addParameter('fignum',[],@isnumeric);
            parser.addParameter('histBins',[],@isnumeric);
            parser.parse(varargin{:});
            p = parser.Results;

            % Plot sampled data to determine how well it fits with the distribution
            if(isempty(p.fignum))
                fig_h = figure();
            else
                fig_h = figure(p.fignum);
            end
            
            if(~isempty(p.plotName))
                fig_h.Name = p.plotName;
            end
            
            clf;
            ax_ws = subplot(2,1,1);
            
            if(isempty(p.time))
                p.time = 1:numel(p.dataSamp);
            end
            
            plot(ax_ws,p.time,p.dataSamp,'b-');
            title(ax_ws,sprintf('%s Timeseries',p.title));
            ylabel(ax_ws,p.dataName);
            xlabel(ax_ws,'Time');
            
            ax_wsh = subplot(2,1,2);
            if ~isempty(p.histBins)
                histogram(ax_wsh,p.dataSamp,p.histBins,'Normalization','probability');
            else
                histogram(ax_wsh,p.dataSamp,'Normalization','probability');
            end
            
            title(ax_wsh,sprintf('%s Distribution',p.title));
            ylabel(ax_wsh,'Probability');
            xlabel(ax_wsh,p.dataName);
        end
    end
    
end