HTM = {
	region = {
			activeColumns = {},
			input = {},
			x,
			y,
			cellsPerColumn,
			desiredLocalActivity,
			permInc,
			permDec,
			minThreshold,
			inhibitionRadius,
			learningOn,
			minOverlap,
			connectedPerm,
			initialPerm,
			minThreshold,
			newSynapseCount,
			overlapRetention,
			activeDutyRetention,
			initialSynapses,
			activationThreshold,
			columns = {
					overlap,
					boost,
					activeDutyCycle,
					overlapDutyCycle,
					active = {_now, _then},
					neighbors = {},
					synapses = {
								index,
								permanence
								},
					cells = {
							activeState = {_now, _then},
							predictiveState = {_now, _then},
							learnState = {_now, _then},
							learningRadius,
							segments = {
										isSequence,
										synapses = {
													column,
													cell,
													permanence
													}
										}
							}
					}
			}
}

function HTM.region:new(_x, _y, _initialRadius, _numSynapses, config)
  config = config or {}
  setmetatable(config, self)
  self.__index = self
  self.x = _x
  self.y = _y
  self.inhibitionRadius = _initialRadius
  self.initialSynapses = _numSynapses
  return config
end

function HTM.region:getNeighbors(_sourceIndex, _radius)
	local _neighbors = {}
	for i = 0, self.x * self.y do
		--loop over every element in the region, return the x,y coords of every cell within _radius of the source
		if i ~= _sourceIndex-1 then
			local _neighborIndex = {}
			local _currX = i % self.y
			local _currY = (i - (i % self.y)) / self.y
				if (math.abs(_currY - (_sourceIndex-1-((_sourceIndex-1) % self.y)) / self.y) <= _radius) and (math.abs(_currX - (_sourceIndex-1) % self.y)) <= _radius then
					table.insert(_neighbors, {_currX, _currY})
				end
		end
	end
	return _neighbors
end

function HTM.region:initialize()
--set the list of neighbors for all columns
--self.columns is ordered according to the dimensions of the region
	for i = 1, self.x * self.y do
		self.columns[i].neighbors = self.getNeighbors(i, self.inhibitionRadius)
		for j = 1, self.initialSynapses do
			local _newIndex = math.random(1, (self.x * self.y))
			local _newPerm = math.random()
			_newPerm = (_newPerm * (.2 * math.random(-1,1))) + self.connectedPerm
			columns[i].synapses[j].index = _newIndex
			columns[i].synapses[j].connectedPerm = _newPerm
		end
	end
end

function HTM.region:spatialPooling()
	--Phase One-------------------------------------------------------------------------------------------------------
	for i = 1, #self.columns do
	--1. for c in columns
	--2.
		self.columns[i].overlap = 0
	--3. overlap(c) = 0
			for j = 1, #self.columns[i].synapses do
				if self.columns[i].synapses[j].permanence > self.connectedPerm then
	--4. for s in connectedSynapses(c)
	--if self.columns[i].synapses[j] permanence is greater than connectedPerm, the synapse is connected.
					local _ti = self.columns[i].synapses[j].index
					if self.inputs[_ti] == 1 then
						self.columns[i].overlap = self.columns[i].overlap + 1
					end
	--5. overlap(c) = overlap(c) + input(t, s.sourceInput)
	--6.
				end
			end
		if self.columns[i].overlap < self.minOverlap then
	--7. if overlap(c) < minOverlap then
			self.columns[i].overlap = 0
	--8. overlap(c) = 0
				else
	--9. else
			self.columns[i].overlap = self.columns[i].overlap * self.columns[i].boost
	--10. overlap(c) = overlap(c) * boost(c)
		end
	end
	--Phase One complete----------------------------------------------------------------------------------------------

	--Phase Two-------------------------------------------------------------------------------------------------------
	for i = 1, #self.columns do
	--11. for c in columns
	--12.
		local _overlaps = {}
		for j = 1, #self.columns[i].neighbors do
			local _index = self.columns[i].neighbors[j]
			table.insert(_overlaps, self.columns[_index].overlap)
		end
		table.sort(_overlaps)
		local _minLocalActivity = _overlaps[self.desiredLocalActivity]
	--13. minLocalActivity = kthScore(neighbors(c), desiredLocalActivity)
	--14.
		if self.columns[i].overlap > 0 and self.columns[i].overlap >= _minLocalActivity then
			table.insert(self.activeColumns, i)
		end
	--15. if overlap(c) > 0 and overlap(c) >= minLocalActivity then
	--16. activeColumns(t).append(c)
	end
	--Phase Two complete.---------------------------------------------------------------------------------------------

	--Phase Three-----------------------------------------------------------------------------------------------------
		for i = 1, #self.activeColumns do
	--18. for c in activeColumns(t)
	--19.
			local activeIndex = self.activeColumns[i]
			for j = 1, #self.columns[activeIndex].synapses do
	--20. for s in potentialSynapses(c)
			local _synapseActiveIndex = self.columns[activeIndex].synapses[j].index
				if self.inputs[_synapseActiveIndex] == 1 then
	--21. if active(s) then
					self.columns[i].synapses[j].permanence = self.columns[i].synapses[j].permanence + self.permInc
	--22. s.permanence += permanenceInc
					if self.columns[i].synapses[j].permanence > 1 then
						self.columns[i].synapses[j].permanence = 1
					end
	--23. s.permanence = min(1.0, s.permanence)
				else
	--24. else
					self.columns[i].synapses[j].permanence = self.columns[i].synapses[j].permanence - self.permDec
	--25. s.permanence -= permanenceDec
					if self.columns[i].synapses[j].permanence < 0 then
						self.columns[i].synapses[j].permanence = 0
					end
	--26. s.permanence = max(0.0, s.permanence)
	--27.
				end
			end
		end
		for i = 1, #self.columns do
	--28. for c in columns:
	--29.
		local _maxDutyCycle = 0
			for j = 1, #self.columns[i].neighbors do
			local _neighborIndex = self.columns[j]
				if self.columns[_neighborIndex].activeDutyCycle > _maxDutyCycle then
					_maxDutyCycle = self.columns[_neighborIndex].activeDutyCycle
				end
			end
		local _minDutyCycle = 0.01 * _maxDutyCycle
	--30. minDutyCycle(c) = 0.01 * maxDutyCycle(neighbors(c))
		self.columns[i].activeDutyCycle = ((self.columns[i].activeDutyCycle * self.activeDutyRetention) + self.inputs[i]) / (self.activeDutyRetention + 1)
	--31. activeDutyCycle(c) = updateActiveDutyCycle(c)
			if self.columns[i].activeDutyCycle > _minDutyCycle then
				self.columns[i].boost = 1
			else
				self.columns[i].boost = _minDutyCycle / self.columns[i].activeDutyCycle
			end
	--32. boost(c) = boostFunction(activeDutyCycle(c), minDutyCycle(c))
	--33.
			if self.columns[i].overlapDutyCycle > minOverlap then
				self.columns[i].overlapDutyCycle = ((self.columns[i].overlapDutyCycle * self.overlapRetention) + 1) / (self.overlapRetention + 1)
			else
				self.columns[i].overlapDutyCycle = (self.columns[i].overlapDutyCycle * self.overlapRetention) / (self.overlapRetention + 1)
			end
	--34. overlapDutyCycle(c) = updateOverlapDutyCycle(c)
		if self.columns[i].overlapDutyCycle < _minDutyCycle then
	--35. if overlapDutyCycle(c) < minDutyCycle(c) then
			for k = 1, #self.columns[i].synapses do
	--36. increasePermanences(c, 0.1*connectedPerm)
	--37.
				self.columns[i].synapses[k].permanence = self.columns[i].synapses[k].permanence + (.1*self.connectedPerm)
			end
		end
	end
		--??? update inhibitionRadius ???
	--Phase Three complete
end

function HTM.region.columns:getPrevActiveSegment(_column, _cell)
	local _segmentActivity = 0
	local _activeSegmentIndex = 0
	for i = 1, #self[_column].cells[_cell].segments do
		for j = 1, #self[_column].cells[_cell].segments[i].synapses do
			local _tempColumnIndex = self[_column].cells[_cell].segments[i].synapses[j].column
			local _tempCellIndex = self[_tempIndex].cells[j].segments[i].synapses[j].cell
			if self.columns[_tempColumnIndex].cells[_tempCellIndex].active._then == 1 then
				_segmentActivity = _segmentActivity + 1
			end
		end
		if _segmentActivity > _maxSegmentActivity then
			if _segmentActivity > self.activationThreshold then
				_maxSegmentActivity = _segmentActivity
				if self[_column].cells[_cell].segments[i].isSequence == 1 then
					_activeSegmentIndex = i
				end
				if self[_column].cells[_cell].segments[i].isSequence == 0 and self[_column].cells[_cell].segments[_activeSegmentIndex].isSequence == 0 then
					_activeSegmentIndex = i
				end
			end
		end
	end
	return(_activeSegmentIndex)
end

function HTM.region.columns:getBestMatchingSegment(_column, _cell)
	local _activeCellIndex = 0
		--get the  segment for cell[a] with the greatest activity, else the most active segment
			local _activeSegmentIndex = 0
			for b = 1, #self[_column].cells[_cell].segments do
				local _segmentActivity = 0
				--iterate over each synapse in the segment and increment _segmentActivity for each active synapse
				for c = 1, #self[_column].cells[_cell].segments[b].synapses do
					local _tempColumnIndex = self[_column].cells[_cell].segments[b].synapses[c].column
					local _tempCellIndex = self[_column].cells[_cell].segments[b].synapses[c].cell
					if self[_tempColumnIndex].cells[_tempCellIndex].active._then == 1 then
						_segmentActivity = _segmentActivity + 1
					end
				end
				if _segmentActivity > _maxSegmentActivity then
					if _segmentActivity > self.minThreshold then
						_maxSegmentActivity = _segmentActivity
						_activeSegmentIndex = b
					end
				end
			end
	return({ _activeSegmentIndex})
end

function HTM.region.columns:getBestMatchingCell(_column)
	local _activeCellIndex = 0
	local _fewestSegments = #self[column].cells[1].segments
	local _fewestSegmentsCell = 0
		for a = 1, #self[_column].cells do
		--get the  segment for cell[a] with the greatest activity, else the most active segment
			local _activeSegmentIndex = 0
			if #self[_column].cells[a].segments < _fewestSegments then
				_fewestSegments = #self[_column].cells[a].segments
				_fewestSegmentsCell = a
			end
			for b = 1, #self[_column].cells[a].segments do
				local _segmentActivity = 0
				--iterate over each synapse in the segment and increment _segmentActivity for each active synapse
				for c = 1, #self[_column].cells[a].segments[b].synapses do
					local _tempColumnIndex = self[_column].cells[a].segments[b].synapses[c].column
					local _tempCellIndex = self[_column].cells[a].segments[b].synapses[c].cell
					if self[_tempColumnIndex].cells[_tempCellIndex].active._then == 1 then
						_segmentActivity = _segmentActivity + 1
					end
				end
				if _segmentActivity > _maxSegmentActivity then
					if _segmentActivity > self.minThreshold then
						_maxSegmentActivity = _segmentActivity
						_activeCellIndex = a
						_activeSegmentIndex = b
					end
				end
			end
		end
		if _activeCellIndex == 0 then
			_activeCellIndex = _fewestSegmentsCell
		end
	return({_activeCellIndex, _activeSegmentIndex})
end

function HTM.region:temporalPooling()
	local sequenceTempUpdateList = {}
	local activeTempUpdateList = {}
	local predictiveTempUpdateList = {}
	local segmentUpdateList = {
								column = {},
								cell = {},
								segment = {
											isSequence,
											synapses = {
														column,
														cell,
														permanence
														}
											}
								}
--Phase One--------------------------------------------------------------------------------------------------------
		for i = 1, #self.activeColumns do
	--18. for c in activeColumns(t)
	--19.
			local _tempIndex = self.activeColumns[i]
			local _bottomUpPredicted = 0
	--20. buPredicted = false
			local _lcChosen = 0
	--21. lcChosen = false
			local _maxSegmentActivity = 0
			for j = 1, #self.columns[_tempIndex].cells do
	--22. for i = 0 to cellsPerColumn - 1
				if self.columns[_tempIndex].cells[j].predictiveState._then == 1 then
	--23. if predictiveState(c, i, t-1) == true then
				local _tempSegmentIndex = self.columns.getPrevActiveSegment(_tempIndex, j)
	--24. s = getActiveSegment(c, i, t-1, activeState)
					if _tempSegmentIndex ~= 0 then
						if self.columns[_tempIndex].cells[j].segments[_tempSegmentIndex].isSequence == 1 then
	--25. if s.sequenceSegment == true then
							_bottomUpPredicted = 1
	--26. buPredicted = true
							self.columns[_tempIndex].cells[j].activeState.now = 1
	--27. activeState(c, i, t) = 1
							if self.columns[_tempIndex].cells[j].segments[_activeSegmentIndex].learnState._then == 1 then
	--28. if segmentActive(s, t-1, learnState) then
								_lcChosen = 1
	--29. lcChosen = true
								self.columns[_tempIndex].cells[j].learnState.now = 1
	--30. learnState(c, i, t) = 1
	--31.
							end
						end
					end
				end
			end
			if _bottomUpPredicted == 0 then
	--32. if buPredicted == false then
				for j2 = 1, #self.columns[_tempIndex].cells do
	--33. for i = 0 to cellsPerColumn - 1
					self.columns[_tempIndex].cells[j2].active.now = 1
	--34. activeState(c, i, t) = 1
	--35.
				end
			end
			if _lcChosen == 0 then
	--36. if lcChosen == false then
				local _bestMatchingCell, _bestMatchingSegment = self.columns.getBestMatchingCell(_tempIndex)[1], self.columns.getBestMatchingCell(_tempIndex)[2]
	--37. I,s = getBestMatchingCell(c, t-1)
				self.columns[_tempIndex].cells[j].learnState.now = 1
	--38. learnState(c, i, t) = 1
				local tempSynapseList = {_column, _cell}
					for d = 1, #self.columns[_tempIndex].cells[_bestMatchingCell].segments[_bestMatchingSegment].synapses do
						local _columnIndex = self.columns[_tempIndex].cells[_bestMatchingCell].segments[_bestMatchingSegment].synapses.column
						local _cellIndex = self.columns[_tempIndex].cells[_bestMatchingCell].segments[_bestMatchingSegment].synapses.cell
						local _cellPerm = self.columns[_tempIndex].cells[_bestMatchingCell].segments[_bestMatchingSegment].synapses.permanence
						if self.columns[_columnIndex].cells[_cellIndex].activeState._then == 1 then
							table.insert(tempSynapseList, {_column = _columnIndex, _cell = _cellIndex, _perm = _cellPerm})
						end
					end
	--39. sUpdate = getSegmentActiveSynapses (c, i, s, t-1, true)
				_tempSegment = {}
				for e = 1, #tempSynapseList do
					table.insert(_tempSegment, {isSequence = 1, column = tempSynapseList[e]._column, cell = tempSynapseList[e]._cell, permanence = tempSynapseList[e]._perm})
	--40. sUpdate.sequenceSegment = true
				end
				table.insert(segmentUpdateList, {column = _tempIndex, cell = j, _tempSegment})
			end
		end
--Phase One end---------------------------------------------------------------------------------------------------------

--Phase Two-------------------------------------------------------------------------------------------------------------
	for i = 1, #self.columns do
		for j = 1, #self.columns[i].cells do
	--42. for c, i in cells
			for k = 1, #self.columns[i].cells[j].segments do
	--43. for s in segments(c, i)
				local segmentActivity = 0
				for l = 1, #self.columns[i].cells[j].segments[k].synapses do
					if self.columns[i].cells[j].activeState.now == 1 then
						segmentActivity = segmentActivity + 1
					end
				end
				if segmentActivity > self.activationThreshold then
	--44. if segmentActive(s, t, activeState) then
					self.columns[i].cells[j].predictiveState.now = 1
	--45. predictiveState(c, i, t) = 1
	--46.
					local tempSynapseList = {_column, _cell, _perm}
					for d = 1, #self.columns[i].cells[j].segments[k].synapses do
						local _columnIndex = sself.columns[i].cells[j].segments[k].synapses.column
						local _cellIndex = self.columns[i].cells[j].segments[k].synapses.cell
						local _cellPerm = self.columns[i].cells[j].segments[k].synapses.permanence
	--47. activeUpdate = getSegmentActiveSynapses (c, i, s, t, false)
						--if synapse to _columnIndex,_cellIndex is active then add it to the potential segment
						if self.columns[_columnIndex].cells[_cellIndex].activeState.now == 1 then
							table.insert(tempSynapseList, {_column = _columnIndex, _cell = _cellIndex, _perm = _cellPerm})
	--48. segmentUpdateList.add(activeUpdate)
	--49.
						end
					end
					local _tempSegment = {}
					for e = 1, #tempSynapseList do
						table.insert(_tempSegment, {isSequence = 0, column = tempSynapseList[e]._column, cell = tempSynapseList[e]._cell, permanence = tempSynapseList[e]._perm})
					end
					table.insert(segmentUpdateList, {column = _tempIndex, cell = j, _tempSegment})
					_bms = self.columnsgetBestMatchingSegment(i,j)
	--50. predSegment = getBestMatchingSegment(c, i, t-1)
					local tempSynapseList2 = {_column, _cell, _perm}
					for d = 1, #self.columns[i].cells[j].segments[_bms].synapses do
						local _columnIndex = sself.columns[i].cells[j].segments[_bms].synapses.column
						local _cellIndex = self.columns[i].cells[j].segments[_bms].synapses.cell
						local _cellPerm = self.columns[i].cells[j].segments[_bms].synapses.permanence
						if self.columns[_columnIndex].cells[_cellIndex].activeState._then == 1 then
							table.insert(tempSynapseList2, {_column = _columnIndex, _cell = _cellIndex, _perm = _cellPerm})
						end
	--51. predUpdate = getSegmentActiveSynapses(c, i, predSegment, t-1, true)
	--52.
					end
					local _tempSegment2 = {}
					for e = 1, #tempSynapseList2 do
						table.insert(_tempSegment2, {isSequence = 0, column = tempSynapseList[e]._column, cell = tempSynapseList[e]._cell, permanence = tempSynapseList[e]._perm})
					end
					table.insert(segmentUpdateList, {column = _tempIndex, cell = j, segment = _tempSegment2})
	--53. segmentUpdateList.add(predUpdate)
				end
			end
		end
	end
--Phase Two end---------------------------------------------------------------------------------------------------------

--Phase Three-----------------------------------------------------------------------------------------------------------
	for i = 1, #self.columns do
		for j = 1, #self.columns[i].cells do
	--54. for c, i in cells
			if self.columns[i].cells[j].learnState._now == 1 then
	--55. if learnState(s, i, t) == 1 then
				for k = 1, #segmentUpdateList do
				local _updateColumn = segmentUpdateList[k].column
				local _updateCell = segmentUpdateList[k].cell
					if _updateColumn == i and _updateCell == j then
						for l = 1, #segmentUpdateList[k].segment do
							segmentUpdateList[k].segment[l].permanence = segmentUpdateList[k].segment[l].permanence + self.permanenceInc
	--56. adaptSegments (segmentUpdateList(c, i), true)
						end
					end
				end
	--57. segmentUpdateList(c, i).delete()
			elseif self.columns[i].cells[j].predictiveState._now == 0 and self.columns[i].cells[j].predictiveState._then == 1 then
	--58. else if predictiveState(c, i, t) == 0 and predictiveState(c, i, t-1)==1 then
				for l = 1, #segmentUpdateList do
				local _updateColumn = segmentUpdateList[l].column
				local _updateCell = segmentUpdateList[l].cell
					if _updateColumn == i and _updateCell == j then
						for m = 1, #segmentUpdateList[l].segment do
							segmentUpdateList[l].segment[m].permanence = segmentUpdateList[l].segment[m].permanence - self.permanenceDec
	--59. adaptSegments (segmentUpdateList(c,i), false)
						end
					end
				end
			end
		end
	end
	--60. segmentUpdateList(c, i).delete()
	--61.
--Phase Three end-------------------------------------------------------------------------------------------------------
end

--
