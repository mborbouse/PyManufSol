from ..GeneralManufSol import *

class GeneralManufSolContainerSinglePhase:
    """ TODO """
    def __init__(self, manuf_sol_list: list[list[GeneralManufSol]], dim: int, spatial_sym_var: list, var_set: VarSetTags = DefaultVarSetTags, temporal_sym_var: list = []):
        """ TODO """
        self.manuf_sol_list_nested = manuf_sol_list
        self.manuf_sol_list = flattenNestedList(manuf_sol_list)
        self.nb_fields = len(self.manuf_sol_list)
        self.dim = dim
        self.var_set = var_set
        self.symUnknownsParams = [list(i.getSymUnknownParam()) for i in self.manuf_sol_list]
        self.symUnknownsParams = flattenNestedList(self.symUnknownsParams)
        sym_var = self.manuf_sol_list[0].getSymIndependentVar()
        for i in range(1,self.nb_fields):
            if self.manuf_sol_list[i].getSymIndependentVar() != sym_var:
                raise ValueError("The independent spatial and temporal coordinates must be identical between the different manuf sols provided.")
        self.spatial_sym_var = spatial_sym_var
        self.temporal_sym_var = temporal_sym_var

    def getNumberFields(self) -> int:
        return self.nb_fields

    def getDim(self) -> int:
        return self.dim

    def getSpatialSymVar(self):
        return self.spatial_sym_var

    def getTemporalSymVar(self):
        return self.temporal_sym_var

    def getVarSetChoice(self) -> VarSetTags:
        return self.var_set

    def getSymUnkownsParams(self) -> list:
        return self.symUnknownsParams

    def getSolsTagsFromVarSetChoice(self, var_set_choice: VarSetTags) -> list[SolutionTags]:
        del var_set_choice
        list_tag = [i for i in DefaultSolutionTags]
        return list_tag[0:self.nb_fields]

    def getSolsTags(self) -> list[SolutionTags]:
        return self.getSolsTagsFromVarSetChoice(self.var_set)
        # """ TODO """
        # tags = [" "] * self.nb_fields
        # for i in range(tags):
        #     tags[i] = str(i)
        # return tags

    def getDicManufSolPerTag(self) -> dict:
        manuf_sol_list_sym = [i.getSymSol() for i in self.manuf_sol_list]
        return assembleManufSolListIntoTagDic(manuf_sol_list_sym, self.getSolsTags())

    def getManufSolFromList(self, ind: int) -> GeneralManufSol:
        return self.manuf_sol_list[ind]

    def changeManufSolCoordSystem(self, currentCoordSystem: CoordinatesSystemType, newSpatialVar: list[sp.Expr], newCoordSystem: CoordinatesSystemType, translation: Optional[list[float]] = None, rotation: Optional[list[list[int,float]]] = None):
        for i in range(len(self.manuf_sol_list_nested)):
            current_spatial_sym_var = self.spatial_sym_var
            new_sol = [j.getSymSol() for j in self.manuf_sol_list_nested[i]]
            if len(self.manuf_sol_list_nested[i]) == 1:
                if currentCoordSystem != newCoordSystem:
                    new_sol[0] = expressScalarInNewCoordSystem(new_sol[0], current_spatial_sym_var, currentCoordSystem, newSpatialVar, newCoordSystem)
                    current_spatial_sym_var = newSpatialVar
                if isinstance(translation, list):
                    new_sol[0] = translateScalarInNewCoordSystem(new_sol[0], current_spatial_sym_var, newSpatialVar, translation)
                    current_spatial_sym_var = newSpatialVar
                if isinstance(rotation, list):
                    new_sol[0] = rotateScalarInNewCoordSystem(new_sol[0], current_spatial_sym_var, newSpatialVar, rotation)
                    current_spatial_sym_var = newSpatialVar    
                self.manuf_sol_list_nested[i][0].setSymSol(new_sol[0])
            else:
                if currentCoordSystem != newCoordSystem:
                    new_sol = expressVectorInNewCoordSystem(new_sol, current_spatial_sym_var, currentCoordSystem, newSpatialVar, newCoordSystem)
                    current_spatial_sym_var = newSpatialVar
                if isinstance(translation, list):
                    new_sol = translateVectorInNewCoordSystem(new_sol, current_spatial_sym_var, newSpatialVar, translation)
                    current_spatial_sym_var = newSpatialVar
                if isinstance(rotation, list):
                    new_sol = rotateVectorInNewCoordSystem(new_sol, current_spatial_sym_var, newSpatialVar, rotation)
                    current_spatial_sym_var = newSpatialVar  
                for j in range(len(self.manuf_sol_list_nested[i])):
                    self.manuf_sol_list_nested[i][j].setSymSol(new_sol[j])

        self.manuf_sol_list = flattenNestedList(self.manuf_sol_list_nested)
        self.spatial_sym_var = newSpatialVar

    # def assembleManufSolListIntoDic(self, manuf_sol_list: list, tags: list[str]) -> dict:
    #     """ TODO """
    #     manuf_sol_dic_per_tag = dict()
    #     tags_for_this_phase = self.getSolsTags()
    #     for tag, isol in zip(tags_for_this_phase,range(self.nb_fields)):
    #         if tag in manuf_sol_dic_per_tag:
    #             manuf_sol_dic_per_tag[tag].append(self.getManufSolFromList(isol).getSymSol())
    #         else:
    #             manuf_sol_dic_per_tag[tag] = self.getManufSolFromList(isol).getSymSol()
    #     return manuf_sol_dic_per_tag

class GeneralManufSolContainerMultiplePhases:
    """ TODO """
    def __init__(self, manuf_sol_lists: list[GeneralManufSolContainerSinglePhase], phase_index_list: list[int]):
        """ TODO """
        if(len(set(phase_index_list)) != len(phase_index_list)):
            raise ValueError("Indices of the phases must be unique.")
        self.nb_phases = len(manuf_sol_lists)
        self.phase_index_list = phase_index_list
        if self.nb_phases != len(phase_index_list):
            raise ValueError("Number of lists of manuf sols must be equal to size of the phase indices list.")
        self.nb_sols_per_phase = dict()
        for sol_phase, iphase in zip(manuf_sol_lists, phase_index_list):
            self.nb_sols_per_phase[iphase] = sol_phase.getNumberFields()
        self.manuf_sol_dic_per_phase = self.assembleManufSolListsIntoDic(manuf_sol_lists)

    def assembleManufSolListsIntoDic(self, manuf_sol_lists: list[GeneralManufSolContainerSinglePhase]) -> dict:
        """ TODO """
        manuf_sol_dic_per_phase = dict()
        for sol_phase, iphase in zip(manuf_sol_lists, self.phase_index_list):
            manuf_sol_dic_per_tag = sol_phase.getDicManufSolPerTag()
            # tags_for_this_phase = sol_phase.getSolsTags()
            # for tag, isol in zip(tags_for_this_phase,range(self.nb_sols_per_phase[iphase])):
            #     if tag in manuf_sol_dic_per_tag:
            #         manuf_sol_dic_per_tag[tag].append(sol_phase.getManufSolFromList(isol).getSymSol())
            #     else:
            #         manuf_sol_dic_per_tag[tag] = sol_phase.getManufSolFromList(isol).getSymSol()
            manuf_sol_dic_per_phase[iphase] = manuf_sol_dic_per_tag
        return manuf_sol_dic_per_phase

    def getDicManufSolPerPhaseAndTag(self):
        return self.manuf_sol_dic_per_phase

    # def __init__(self, manuf_sol_lists: list[GeneralManufSol]):
    #     """ TODO """
    #     self.manuf_sol_lists = manuf_sol_lists
    #     self.nb_phases = len(manuf_sol_lists)
    #     self.nb_sols_per_phase = [0] * self.nb_phases
    #     for sol_phase, iphase in zip(manuf_sol_lists, range(self.nb_phases)):
    #         self.nb_sols_per_phase[iphase] = len(sol_phase)

    # def getSolsTagsForThisPhase(self, iphase: int) -> list[str]:
    #     """ TODO """
    #     tags = [" "] * self.nb_sols_per_phase[iphase]
    #     for i in range(self.nb_sols_per_phase[iphase]):
    #         tags[i] = str(i)
    #     return tags

    # def assembleManufSolListsIntoDic(self) -> dict:
    #     """ TODO """
    #     manuf_sol_dic_per_phase = dict()
    #     for sol_phase, iphase in zip(self.manuf_sol_lists, range(self.nb_phases)):
    #         manuf_sol_dic_per_tag = dict()
    #         tags_for_this_phase = self.getSolsTagsForThisPhase(iphase)
    #         for tag, isol in zip(tags_for_this_phase,range(self.nb_sols_per_phase[iphase])):
    #             if tag in manuf_sol_dic_per_tag:
    #                 manuf_sol_dic_per_tag[tag].append(sol_phase[isol].getSymSol())
    #             else:
    #                 manuf_sol_dic_per_tag[tag] = sol_phase[isol].getSymSol()
    #         manuf_sol_dic_per_phase[iphase] = manuf_sol_dic_per_tag
    #     return manuf_sol_dic_per_phase