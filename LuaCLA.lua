local region = {}

function region:new(x, y, initialRadius, numSynapses, cellsPerColumn, activeDutyRetention, minOverlap, overlapRetention, initialPerm, connectedPerm)
	local config = {}
	self.x = x
	self.y = y
	self.cellsPerColumn = cellsPerColumn
	self.connectedPerm = 3
	self.inhibitionRadius = initialRadius
	self.initialSynapses = numSynapses
	self.activeDutyRetention = activeDutyRetention
	self.minOverlap = minOverlap
	self.overlapRetention = overlapRetention
	self.input = {}
	self.columns = {}
	self.initialPerm = initialPerm
	self.connectedPerm = connectedPerm
	self.length = x*y
	print(self.length)

	for i = 1, self.length do
		self.columns[i] = {	overlap = 0,
							boost = 0,
							active = 0,
							minDutyCycle = 0,
							activeDutyCycle = 0,
							overlapDutyCycle = 0,
							isActive = 0,
							wasActive =0,
							neighbors = {},
							synapses = {},
							cells = {}
							}
	end

	for _, column in ipairs(self.columns) do
		for _, cell in ipairs(column.cells) do
			cell = {	isActive = 0,
						wasActive = 0,
						isPredicting = 0,
						wasPredicting = 0,
						isLearning = 0,
						wasLearning = 0,
						segments = {
										sequenceSegment  = 0,
										activationThreshold,
										synapses = {},
										segmentUpdateList = {}
									}
					}
		end

		for i = 1, numSynapses do
			table.insert(column.synapses, {permanence = self.initialPerm, column = math.random(self.length), cell = math.random(cellsPerColumn)})
		end
	end

	setmetatable(config, self)
	self.__index = self
return config
end

function region:getNeighbors(sourceIndex, radius)
	local neighbors = {}
	for i = 0, (self.x * self.y) do
		--loop over every element in the region, return the x,y coords of every cell within _radius of the source
		if i ~= sourceIndex-1 then
			local neighborIndex = {}
			local currX = i % self.y
			local currY = (i - (i % self.y)) / self.y
				if (math.abs(currY - (sourceIndex-1-((sourceIndex-1) % self.y)) / self.y) <= radius) and (math.abs(currX - (sourceIndex-1) % self.y)) <= radius then
					table.insert(neighbors, {currX, currY})
				end
		end
	end
	return neighbors
end

function activeColumns(t)
	local function activeColumns_it(t, i)
		i = i + 1
		while t[i] and t[i].active == 0 do i = i + 1 end
		if t[i] then return t[i] end
	end
	return activeColumns_it, t, 0
end

function getConnectedSynapses(t, connectedPerm)
	local synapseTable = {}
	for i = 1, #t.synapses do
		if t.synapses[i].permanence > connectedPerm then
			table.insert(synapseTable, i)
		end
	end
	return #synapseTable
end

function region:spatialPooling()
	local columns = self.columns
	local numColumns = #columns
	local connectedPerm = self.connectedPerm
	local minOverlap = self.minOverlap
	local desiredLocalActivity = self.desiredLocalActivity
	local permanenceInc = self.permanenceInc
	local permanenceDec = self.permanenceDec
	local activeDutyRetention = self.activeDutyRetention
	local inputs = self.inputs
	local overlapRetention = self.overlapRetention
	local inhibitionRadius = self.inhibitionRadius

	for _, column in ipairs(columns) do
		local numConnected = 0
		for _, synapse in ipairs(column.synapses) do
			if synapse.permanence > connectedPerm then
				numConnected = numConnected + 1
			end
		end

		column.overlap = numConnected
		if column.overlap < minOverlap then
				column.overlap = 0
		else
				column.overlap = column.overlap * column.boost
		end
	end
	for _,column in ipairs(columns) do
		local overlaps = {}
		for _,neighbor in ipairs(column.neighbors) do
			table.insert(overlaps, neighbor.index)
		end
		table.sort(overlaps)
		local minLocalActivity = overlaps[desiredLocalActivity]
		if column.overlap > 0 and column.overlap >= minLocalActivity then
			column.isActive = 1
		end
	end

	print(activeColumns)
	for column in activeColumns(columns) do
		for _, synapse in ipairs(column.synapses) do
			local cellIndex = synapse.cell
			local columnIndex = synapse.column
			if self.columns[columnIndex].cell[cellIndex].isActive == 1 then
				synapse.permanence = synapse.permanence + permanenceInc
				if synapse.permanence > 1 then synapse.permanence = 1 end
			else
				synapse.permanence = synapse.permanence - permanenceDec
				if synapse.permanence < 0 then synapse.permanence = 0 end
			end
		end
	end

	for columnIndex, column in ipairs(columns) do
		local maxDutyCycle = 0
		for _, neighbor in ipairs(column.neighbors) do
			if neighbor.activeDutyCycle > maxDutyCycle then
				maxDutyCycle = neighbor.activeDutyCycle
			end
		end
		column.minDutyCycle = .01 * maxDutyCycle
		column.activeDutyCycle = ((column.activeDutyCycle * activeDutyRetention) + self.input[columnIndex]) / (activeDutyRetention + 1)
		minDutyCycle = column.minDutyCycle
		if column.activeDutyCycle > minDutyCycle then
			column.boost = 1
		else
			column.boost = minDutyCycle / column.activeDutyCycle
		end
		if column.overlapDutyCycle > minOverlap then
			column.overlapDutyCycle = ((column.overlapDutyCycle * overlapRetention) + 1) / (overlapRetention + 1)
		else
			column.overlapDutyCycle = (column.overlapDutyCycle * overlapRetention) / (overlapRetention + 1)
		end
		if column.overlapDutyCycle < column.minDutyCycle then
			for _, synapse in ipairs(column.synapses) do
				synapse.permanence = synapse.permanence + (.1*connectedPerm)
			end
		end
		--inhibitionRadius = averageReceptiveFieldSize()
	end
end

function region:temporalPooling()
	local activationThreshold = self.activationThreshold
	local columns = self.columns
	for column in activeColumns(columns) do
		local buPredicted = 0
		local lcChosen = 0
		for _,cell in ipairs(column.cells) do
			if cell.wasPredicting == 1 then
				local seqSegMax = 0
				local seqSegIndex = 0
				local segMax = 0
				local segIndex = 0
				local segIsLearning = 0
				for _, segment in ipairs(cell.segments) do
					local tempActivity = 0
					local tempLearning = 0
					for _, synapse in ipairs(segment.synapses) do
						local colIndex = synapse.column
						local celIndex = synapse.cell
						if columns[colIndex].cells[celIndex].wasActive then tempActivity = tempActivity + 1 end
						if columns[colIndex].cells[celIndex].wasLearning then tempLearning = tempLearning + 1 end
					end
					if segment.isSequence then
						seqSegMax = tempActivity
						seqSegIndex = k
						segMax = tempActivity
						segIndex = k
						if tempLearning > segment.activationThreshold then
							segIsLearning = 1
						else
							segIsLearning = 0
						end
					else if seqSegIndex == 0 and tempActivity > segMax then
						segMax = tempActivity
						segIndex = k
						if tempLearning > segment.activationThreshold then
							segIsLearning = 1
						else
							segIsLearning = 0
						end
					end
				end
				if seqSegIndex == 1 then
					buPredicted = 1
					cell.isActive = 1
					if segIsLearning == 1 then
						lcChosen = 1
						cell.isLearning = 1
					end
				end
			end
		end
		if buPredicted == 0 then
			for _, cell in ipairs(column.cells) do
				cell.isActive = 1
			end
		end
		if lcChosen == 0 then
			local bestCell = 0
			local segMax = 0
			local segmentIndex = 0
			for cellIndex, cell in ipairs(column.cells) do
				for segmentIndex_inner, segment in ipairs(cell.segments) do
					local tempActivity = 0
					for _, synapse in ipairs(segment.synapses) do
						if columns[synapse.column].cells[synapse.cell].wasActive and synapse.permanence > 0 then
							tempActivity = tempActivity + 1
						end
					end
					if tempActivity > segMax then
						segMax = tempActivity
						segmentIndex = segmentIndex_inner
						bestCell = cellIndex
					end
				end
			end

			if bestCell == 0 then
				local tempLeastSegs = math.huge
				local tempLeastSegsCell = 0
				for _, cell in ipairs(column.cells) do
					local tempNumSegs = #cell.segments
					if tempNumSegs < tempLeastSegs then
						tempLeastSegs = tempNumSegs
					end
				end
				bestCell = tempLeastSegsCell
			end
			column.cells[bestCell].isLearning = 1
			local sUpdate = {}
			sUpdate.synapses = {}
			sUpdate.column = i
			sUpdate.cell = bestCell
			if segmentIndex ~= 0 then
				for _, synapse in ipairs(cell.segments[segmentIndex].synapses) do
					local tempActivity = 0
						local columnIndex = synapse.column
						local cellIndex = synapse.cell
						if columns[columnIndex].cells[cellIndex].wasActive > 0 then
							table.insert(sUpdate.synapses, synapse)
						end
					end
				end
			end
			sUpdate.sequenceSegment = 1
			table.insert(columns.cells[bestCell].segments[segmentIndex].segmentUpdateList, sUpdate)
		end
	end

	for _, column in ipairs(columns) do
		for _, cell in ipairs(column.cells) do
			for _, segment in ipairs(cell.segments) do
				local tempActivity = 0
				for _, synapse in ipairs(segment.synapses) do
					local columnIndex = synapse.column
					local cellIndex = synapse.cell
					if columns[columnIndex].cells[cellIndex].isActive then
						tempActivity = tempActivity + 1
					end
				end
				local segIsActive = 0
				if tempActivity > activationThreshold then
					segIsActive = 1
				end
				if segIsActive == 1 then
					columns[columnIndex].cells[cellIndex].isPredicting = 1
					local activeUpdate = {}
					sUpdate.synapses = {}
					for _, synapse in ipairs(segment.synapses) do
						local columnIndex = synapse.column
						local cellIndex = synapse.cell
						if columns[columnIndex].cells[cellIndex].isActive then
							table.insert(activeUpdate, segment)
						end
					end
					table.insert(segment.segmentUpdateList, activeUpdate)
					local segMax = 0
					local segIndex = 0
					local tempActivity = 0
					for segmentIndex, segment in ipairs(columns[columnIndex].cells[cellIndex].segments) do
						for _, synapse in ipairs(segment.synapses) do
							local columnIndex_in = synapse.column
							local cellIndex_in = synapse.cell
							if columns[columnIndex_in].cells[cellIndex_in].wasPredicting then
								tempActivity = tempActivity + 1
							end
						end
						if tempActivity > segMax then
							segMax = tempActivity
							segIndex = segmentIndex
						end
					end
					local predSegment = cell.segments[segmentIndex]
					local activeUpdateSynapses = {}
					for _,synapse in ipairs(predSegment) do
						if columns[synapse.column].cells[synapse.cell].isActive then
							table.insert(activeUpdateSynapses, synapse)
						end
					end
					table.insert(segment.segmentUpdateList, activeUpdateSynapses)
				end
			end
		end
	end

	for _, column in ipairs(columns) do
		for _, cell in ipairs(column.cells) do
			for _, segment in ipairs(cell.segments) do
				for _, update in ipairs(segment.segmentUpdateList) do
					if cell.isLearning == 1 then
						for _, synapse in ipairs(segment.synapses) do
							if update.column == synapse.column then
								if update.cell == synapse.cell then
								synapse.permanence = synapse.permanence + permanenceInc
								end
							else
								update.permanence = initialPerm
								table.insert(segment, update)
							end
						end
					elseif cell.wasLearning == 1 and cell.isLearning == 0 then
						for _, synapse in ipairs(segment.synapses) do
							if update.column == synapse.column then
								if update.cell == synapse.cell then
								synapse.permanence = synapse.permanence - permanenceDec
								end
							else
								update.permanence = initialPerm
								table.insert(segment, update)
							end
						end
					end
					update = nil
				end
			end
		end
	end
end

function region:cycleStates()
	for _, column in ipairs(self.columns) do
		column.wasActive = 0
		if column.isActive == 1 then
			column.isActive = 0
			column.wasActive = 1
		end
		for _, cell in ipairs(column.cells) do
			cell.wasActive = 0
			if cell.isActive == 1 then
				cell.isActive = 0
				cell.wasActive = 1
			end
		end
	end
end

--testRegion is a basic region, and is fed a single input state to test if the damn thing even runs.

--region:new(x, y, initialRadius, numSynapses, cellsPerColumn, activeDutyRetention, minOverlap, overlapRetention, initialPerm, connectedPerm, initialSynapses)
testRegion = region:new(10, 10 , 4, 6, 4, 1000, 5, 1000, .4, .5, 12)

testRegion.input = {1,0,0,0,0,0,0,0,0,0,
					1,0,0,0,0,0,0,0,0,0,
					1,0,0,0,0,0,0,0,0,0,
					1,0,0,0,0,0,0,0,0,0,
					1,0,0,0,0,0,0,0,0,0,
					1,0,0,0,0,0,0,0,0,0,
					1,0,0,0,0,0,0,0,0,0,
					1,0,0,0,0,0,0,0,0,0,
					1,0,0,0,0,0,0,0,0,0,
					1,0,0,0,0,0,0,0,0,0}

testRegion:spatialPooling()
testRegion:temporalPooling()
testRegion:cycleStates()
