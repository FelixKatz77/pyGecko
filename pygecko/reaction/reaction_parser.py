import pandas as pd
from datetime import datetime
from pygecko.reaction.well_plate import Well_Plate
from ord_schema.proto import reaction_pb2
from ord_schema.proto import dataset_pb2
from ord_schema.units import UnitResolver
from ord_schema import validations
from ord_schema import message_helpers
unit_resolver = UnitResolver()

class Reaction_Parser:

    '''
    A class wrapping the functionality of the ord-schema library of the open reaction database to create a dataset from
    a combinatorial reaction layout.
    '''

    @classmethod
    def build_dataset(cls, layout: Well_Plate, results_df: pd.DataFrame,
                      path: str|None = None) -> dataset_pb2.Dataset:

        '''
        Returns a ORD dataset created from a combinatorial reaction layout and writes it to a file if a path is given.

        Args:
            layout (Well_Plate): Well_Plate object containing the combinatorial reaction layout.
            results_df (pd.DataFrame): DataFrame containing the results of the reactions.
            path (str|None, optional): Path to write the dataset to. Defaults to None.

        Returns:
            dataset_pb2.Dataset: ORD dataset created from the combinatorial reaction layout.
        '''

        dataset = cls.create_dataset_from_layout(layout, results_df)
        if path:
            message_helpers.write_message(dataset, path)
        return dataset

    @classmethod
    def create_dataset_from_layout(cls, layout:Well_Plate, results_df: pd.DataFrame) -> dataset_pb2.Dataset:

        '''
        Returns a ORD dataset created from a combinatorial reaction layout.

        Args:
            layout (Well_Plate): Well_Plate object containing the combinatorial reaction layout.
            results_df: DataFrame containing the results of the reactions.

        Returns:
            dataset_pb2.Dataset: ORD dataset created from the combinatorial reaction layout.
        '''

        reactions = []
        for x in layout.design.x.values:
            for y in layout.design.y.values:
                _yield = results_df.loc[x, str(y)]
                reaction = reaction_pb2.Reaction()
                reaction.identifiers.add(value=f"{x}{y}", type="CUSTOM", details="Plate position")
                cls.__add_inputs(reaction, layout.meta_data, x, y)
                cls.__add_conditions(reaction, layout.meta_data)
                cls.__add_workups(reaction, layout.meta_data)
                cls.__add_outcomes(reaction, layout.meta_data, _yield, layout.get_product(f"{x}{y}"))
                cls.__add_provenance(reaction, layout.meta_data)
                reactions.append(reaction)
        dataset = dataset_pb2.Dataset(reactions=reactions)
        validations.validate_message(dataset)
        return dataset


    @classmethod
    def __add_inputs(cls, reaction:reaction_pb2.Reaction, metadata:dict, x:str, y:str) -> None:

        '''
        Adds the inputs to a reaction.

        Args:
            reaction (reaction_pb2.Reaction): Reaction to add the inputs to.
            metadata (dict): Metadata of the reaction.
            x (str): x-coordinate of the reaction in the combinatorial reaction layout.
            y (str): y-coordinate of the reaction in the combinatorial reaction layout.
        '''

        for stock in metadata['stock_solutions']:
            if x in stock['wells'] or str(y) in stock['wells'] or stock['wells'] == 'all':
                solute = reaction.inputs[stock["name"]].components.add()
                solvent = reaction.inputs[stock["name"]].components.add()
                inputs = cls.__create_inputs_from_stock(stock)
                solute.CopyFrom(inputs[0])
                solvent.CopyFrom(inputs[1])


    @staticmethod
    def __add_conditions(reaction:reaction_pb2.Reaction, metadata:dict) -> None:

        '''
        Adds the conditions to a reaction.

        Args:
            reaction (reaction_pb2.Reaction): Reaction to add the conditions to.
            metadata (dict): Metadata of the reaction.
        '''

        conditions = metadata["conditions"]
        reaction.setup.vessel.CopyFrom(
            reaction_pb2.Vessel(
                type=conditions["vessel"]["type"],
                material=dict(type=conditions["vessel"]["material"], details=conditions["vessel"]["model"]),
                volume=unit_resolver.resolve(conditions["vessel"]["volume"]),
        ))
        reaction.conditions.temperature.setpoint.CopyFrom(unit_resolver.resolve(conditions['temperature']))

    @staticmethod
    def __add_workups(reaction:reaction_pb2.Reaction, metadata:dict) -> None:

        '''
        Adds the workup steps to a reaction.

        Args:
            reaction (reaction_pb2.Reaction): Reaction to add the workup steps to.
            metadata (dict): Metadata of the reaction.
        '''

        for step in metadata["workup"]:
            workup = reaction.workups.add()
            workup.CopyFrom(
                reaction_pb2.ReactionWorkup(
                    type=step["type"],
                    keep_phase=step["keep_phase"],
                    duration = unit_resolver.resolve(step["duration"]),
                    details=step["details"]
                ))

    @staticmethod
    def __add_outcomes(reaction:reaction_pb2.Reaction, metadata:dict, _yield:float, product_mol:str) -> None:

        '''
        Adds the outcomes to a reaction.

        Args:
            reaction (reaction_pb2.Reaction): Reaction to add the outcomes to.
            metadata (dict): Metadata of the reaction.
            _yield (float): Yield of the reaction.
            product_mol (str): SMILES string of the product of the reaction.
        '''

        conditions = metadata["conditions"]
        outcome = reaction.outcomes.add()
        outcome.reaction_time.CopyFrom(unit_resolver.resolve(conditions["time"]))
        for analysis in metadata["analysis"]:
            if analysis["purpose"] == "QUANTIFICATION":
                quantitative_analysis = analysis
            outcome.analyses[analysis["type"]].CopyFrom(
                reaction_pb2.Analysis(
                    type=analysis["type"],
                    details=analysis["details"],
                    instrument_manufacturer=analysis["instrument_manufacturer"]))
        if _yield:
            product = outcome.products.add()
            product.identifiers.add(value=product_mol, type="SMILES")
            product.is_desired_product = True
            product.reaction_role = reaction_pb2.ReactionRole.PRODUCT
            measurement = product.measurements.add(analysis_key=quantitative_analysis["type"], type="AREA")
            measurement.is_normalized = quantitative_analysis["is_normalized"]  # peak areas are relative
            measurement.percentage.value = _yield  # placeholder

    @staticmethod
    def __add_provenance(reaction:reaction_pb2.Reaction, metadata:dict) -> None:

        '''
        Adds the provenance to a reaction.

        Args:
            reaction (reaction_pb2.Reaction): Reaction to add the provenance to.
            metadata (dict): Metadata of the reaction.
        '''

        provenance = metadata["provenance"]
        reaction.provenance.city = provenance["city"]
        reaction.provenance.doi = provenance["doi"]
        reaction.provenance.record_created.time.value = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
        reaction.provenance.record_created.person.CopyFrom(
            reaction_pb2.Person(
                name=provenance["author"],
                organization=provenance["organization"],
                orcid=provenance["orcid"],
                email=provenance["email"]))

    @staticmethod
    def __create_inputs_from_stock(stock:dict) -> list[reaction_pb2.Compound]:

        '''
        Creates the inputs for a reaction from a stock solution.

        Args:
            stock (dict): Dictionary describing a stock solution to create the inputs from.

        Returns:
            list[reaction_pb2.Compound]: List of compounds to add to the reaction.
        '''

        compound = message_helpers.build_compound(
            smiles=stock['compound'],
            role=stock['role']
        )
        solvent = message_helpers.build_compound(
            name=stock['solvent'],
            role='SOLVENT',
        )
        solvent.amount.volume.CopyFrom(unit_resolver.resolve(stock['volume']))
        input = message_helpers.set_solute_moles(
            solute=compound,
            concentration=stock["concentration"],
            solvents=[solvent]
        )
        return input
