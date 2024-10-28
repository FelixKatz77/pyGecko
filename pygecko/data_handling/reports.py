import io
from functools import partial
from datetime import datetime
from pathlib import Path

import numpy as np
from pygecko.visualization.visuals import Visualization
from reportlab.graphics.shapes import Drawing, Line
from reportlab.lib import utils
from reportlab.lib.colors import Color
from reportlab.lib.pagesizes import A4, portrait
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, Frame, PageTemplate, \
    FrameBreak, NextPageTemplate, PageBreak
from reportlab.lib.units import cm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from indigo import *
from indigo.renderer import IndigoRenderer
from pygecko.reaction import Reaction_Array

indigo = Indigo()
renderer = IndigoRenderer(indigo)
styles=getSampleStyleSheet()
report_style = ParagraphStyle('report_style', alignment=0, allowOrphans=0, allowWidows=1, backColor=None, borderColor=None,
                            borderPadding=0, borderRadius=None, borderWidth=0, firstLineIndent=0, fontName='Helvetica',
                            fontSize=10, spaceAfter=4, spaceBefore=6, textColor=Color(0,0,0,1), leading=12)
report_style_bold = ParagraphStyle('report_style', alignment=0, allowOrphans=0, allowWidows=1, backColor=None, borderColor=None,
                                borderPadding=0, borderRadius=None, borderWidth=0, firstLineIndent=0, fontName='Helvetica-Bold',
                                fontSize=10, spaceAfter=4, spaceBefore=6, textColor=Color(0.1,0.1,0.44,1), leading=12)

class Report:

    def __init__(self, experiment_id:str, layout:Reaction_Array, yield_array:np.ndarray):
        self.id = experiment_id
        self.layout = layout
        self.yield_array = yield_array
        pass

class PDF_Report(Report):

    def __init__(self, filename: str, experiment_id:str, layout:Reaction_Array, yield_array:np.ndarray):
        super().__init__(experiment_id, layout, yield_array)
        self.doc = SimpleDocTemplate(filename, pagesize=A4, rightMargin=1.5*cm, leftMargin=1.5*cm, topMargin=2*cm,
                                    bottomMargin=2*cm, showBoundary=0)
        self.author = self.layout.meta_data['provenance']['author']
        self.analysis = ', '.join(analysis['type'] for analysis in self.layout.meta_data['analysis'])
        self.metadata = self.__collect_meta_data()
        self.analysis_table = self.__create_analysis_table()
        heatmap_path = Path(__file__).resolve().parent.joinpath('tmp/heatmap.png')
        Visualization.visualize_plate(yield_array, results='yield', path=heatmap_path)
        self.heatmap = Image(heatmap_path, 11*cm, 6*cm)
        self.results_table = self.__create_results_table_quantification()
        self.molecules_table_rows, self.molecules_table_columns, self.common_molecules_img = self.__create_molecules_table()
        self.conditions_table = self.create_conditions_table()
        #if self.analysis == 'Quantification':
            #self.product_structure = self.__create_product_image(report_data.analysis_data['Product'][1][0])
        #else:
        self.product_structure = None
        self.__create_pdf_report()


    def __create_pdf_report(self):
        elements = []
        h_padding = 15
        h_w, h_h = self.heatmap.wrap(0,0)
        rt_w, rt_h = self.results_table.wrap(0, 0)
        p_w, p_h = self.metadata.wrap(self.doc.width / 3, self.doc.height)
        cd_w, cd_h = self.conditions_table.wrap(0, 0)
        cm_w, cm_h = self.analysis_table.wrap(0, 0)
        p_h += h_padding
        cm_h += 25 + h_padding
        cd_h += 25 + h_padding
        h_h += 25 + h_padding
        rt_h += 25 + h_padding
        mt_w, mt_h = self.molecules_table_columns.wrap(0,0)
        mt_h += 20 + h_padding
        com_w, com_h = self.common_molecules_img.wrap(0,0)
        com_h += 20 + h_padding

        if self.analysis == 'Filter':
            rt_h = rt_h/2

        header = partial(self.__create_header)

        top_fr_l = Frame(self.doc.leftMargin, self.doc.height + self.doc.bottomMargin - p_h, p_w, p_h, leftPadding=0,
                        bottomPadding=5, rightPadding=0, topPadding=0, showBoundary=0, id='f1')
        top_fr_r = Frame(self.doc.leftMargin + p_w, self.doc.height + self.doc.bottomMargin - p_h, self.doc.width - p_w,
                        p_h, leftPadding=0, bottomPadding=5, rightPadding=0, topPadding=0, showBoundary=0, id='f2')
        conditions_fr = Frame(self.doc.leftMargin, self.doc.height + self.doc.bottomMargin - p_h - cd_h, self.doc.width,
                              cd_h, leftPadding=0, bottomPadding=0, rightPadding=0, topPadding=0, showBoundary=0, id='f3')
        comments_fr = Frame(self.doc.leftMargin,
                            self.doc.height + self.doc.bottomMargin - p_h -  cd_h - cm_h,
                            self.doc.width, cm_h, leftPadding=0, bottomPadding=0, rightPadding=0, topPadding=5,
                            showBoundary=0, id='f4')
        results_fr = Frame(self.doc.leftMargin, self.doc.height + self.doc.bottomMargin - p_h - cd_h - cm_h - rt_h,
                           self.doc.width,
                           rt_h, leftPadding=0, bottomPadding=0, rightPadding=0, topPadding=0,
                           showBoundary=0, id='f6')
        heatmap_fr = Frame(self.doc.leftMargin, self.doc.height + self.doc.bottomMargin - p_h - cd_h - cm_h - rt_h - h_h,
                           self.doc.width, h_h, leftPadding=0, bottomPadding=0, rightPadding=0, topPadding=0,
                           showBoundary=0, id='f4')
        setup_fr = Frame(self.doc.leftMargin, self.doc.height + self.doc.bottomMargin - 25, self.doc.width, 25,
                        leftPadding=0, bottomPadding=0, rightPadding=0, topPadding=0, showBoundary=0, id='f7')
        molecules_table_l_fr = Frame(self.doc.leftMargin, self.doc.height + self.doc.bottomMargin - 25 - mt_h, self.doc.width/2, mt_h,
                                    leftPadding=0, bottomPadding=0, rightPadding=0, topPadding=0, showBoundary=0, id='f8')
        molecules_table_r_fr = Frame(self.doc.leftMargin+self.doc.width/2, self.doc.height + self.doc.bottomMargin - 25 - mt_h, self.doc.width/2, mt_h,
                                    leftPadding=0, bottomPadding=0, rightPadding=0, topPadding=0, showBoundary=0, id='f9')
        common_molecules_fr = Frame(self.doc.leftMargin, self.doc.height + self.doc.bottomMargin - 25 - mt_h - com_h, self.doc.width, com_h, leftPadding=0,
                                    bottomPadding=5, rightPadding=0, topPadding=0, showBoundary=0,
                                    id='f10')
        thin_line = Drawing(0, 0)
        thin_line.add(Line(0, 0, self.doc.width, 0, strokeWidth=0.8, strokeColor=Color(0.1,0.1,0.44,1)))

        thick_line = Drawing(0, 0)
        thick_line.add(Line(0, 0, self.doc.width, 0, strokeWidth=1))


        elements.append(self.metadata)
        elements.append(FrameBreak)
        if self.product_structure:
            elements.append(Paragraph('Product:', report_style_bold))
            elements.append(self.product_structure)
        elements.append(FrameBreak)
        #elements.append(thick_line)
        elements.append(Spacer(width=0, height=5))
        elements.append(Paragraph('Conditions:', report_style_bold))
        elements.append(thin_line)
        elements.append(Spacer(width=0, height=5))
        elements.append(self.conditions_table)
        elements.append(Spacer(width=0, height=5))
        elements.append(FrameBreak)
        elements.append(Paragraph('Analysis:', report_style_bold))
        elements.append(thin_line)
        elements.append(Spacer(width=0, height=5))
        elements.append(self.analysis_table)
        elements.append(Spacer(width=0, height=5))
        elements.append(FrameBreak)
        elements.append(Paragraph('Results:', report_style_bold))
        elements.append(thin_line)
        elements.append(Spacer(width=0, height=5))
        elements.append(self.results_table)
        #elements.append(Spacer(width=0, height=5))
        elements.append(FrameBreak)
        elements.append(Paragraph('Heatmap:', report_style_bold))
        elements.append(thin_line)
        elements.append(Spacer(width=0, height=5))
        elements.append(self.heatmap)
        elements.append(NextPageTemplate('Second'))
        elements.append(PageBreak())
        elements.append(Paragraph('Plate Setup:', report_style_bold))
        elements.append(Spacer(width=0, height=2))
        elements.append(thin_line)
        elements.append(FrameBreak)
        elements.append(Paragraph('Rows:', report_style))
        elements.append(self.molecules_table_rows)
        elements.append(FrameBreak)
        elements.append(Paragraph('Columns:', report_style))
        elements.append(self.molecules_table_columns)
        elements.append(FrameBreak)
        elements.append(thin_line)
        elements.append(Spacer(width=0, height=4))
        elements.append(Paragraph('All Wells:', report_style))
        elements.append(self.common_molecules_img)

        #if self.analysis == 'Quantification':
        self.doc.addPageTemplates([PageTemplate(id='First', frames=[top_fr_l, top_fr_r, conditions_fr, comments_fr, results_fr, heatmap_fr], onPage=header)])
        # else:
        #     self.doc.addPageTemplates([PageTemplate(id='First',
        #                                             frames=[top_fr_l, top_fr_r, conditions_fr, results_txt_fr,
        #                                                     results_fr1, results_fr2, heatmap_fr, comments_fr], onPage=header)])
        self.doc.addPageTemplates([PageTemplate(id='Second', frames=[setup_fr, molecules_table_l_fr, molecules_table_r_fr, common_molecules_fr], onPage=header)])

        self.doc.build(elements)

    def __collect_meta_data(self):
        metadata = f'<b>ID: {self.id}</b>'
        metadata += f'<br />Author: {self.author}'
        metadata += f'<br />Analysis: {self.analysis}'
        p_metadata = Paragraph(metadata, report_style)
        return p_metadata


    # TODO: values in report_df that correspond to unused wells must be Nan

    def __create_results_table_quantification(self):
        alignment = 'CENTRE'
        yield_list = []
        for lst in self.yield_array.tolist():
            row = []
            for element in lst:
                if element == -1:
                    element = 'n.d.'
                row.append(element)
            yield_list.append(row)
        df_list = [list(self.layout.columns.keys())] + yield_list
        for i in range(len(df_list) - 1):
            df_list[i + 1].insert(0, list(self.layout.rows.keys())[i])
        df_list[0].insert(0, '')
        table = Table(df_list, colWidths=1.2 * cm, hAlign=alignment)
        table.setStyle(TableStyle([('ALIGNMENT', (0, 0), (-1, -1), 'CENTRE'),
                                   ('FONTSIZE', (0, 0), (-1, -1), 9), ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                                   ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
                                   ('FONTNAME', (1, 1), (-1, -1), 'Helvetica')]))
        return table

    # def __create_results_table_filtering(self, report_data:Report_Data):
    #     if self.analysis == 'Quantification':
    #         alignment = 'CENTRE'
    #     else:
    #         alignment = 'LEFT'
    #     report_df = report_data.report_df
    #     df_list = [report_df.columns[:, ].values.astype(str).tolist()] + report_df.values.tolist()
    #     for i in range(len(df_list) - 1):
    #         df_list[i + 1].insert(0, report_df.index[i])
    #     df_list[0].insert(0, '')
    #     table = Table(df_list, colWidths=1.2 * cm, hAlign=alignment)
    #     table.setStyle(TableStyle([('ALIGNMENT', (0, 0), (-1, -1), 'CENTRE'),
    #                             ('FONTSIZE', (0, 0), (-1, -1), 9), ('FONTNAME', (0, 0), (-1, 0), 'Times-Bold'),
    #                             ('FONTNAME', (0, 0), (0, -1), 'Times-Bold'),
    #                             ('FONTNAME', (1, 1), (-1, -1), 'Times-Roman')]))
    #     return table

    def __create_header(self, canvas, doc):
        canvas.saveState()

        canvas.setFont('Helvetica', 9)
        (width, height) = portrait(A4)
        timestamp = datetime.now()
        date = timestamp.strftime(('%d.%m.%Y %H:%M'))
        canvas.drawString(1*cm, height - 1*cm , f"pyGECKO Report - {self.id}")
        canvas.drawRightString(width - 1*cm, height - 1*cm, date)

        canvas.restoreState()

    def __create_molecules_table(self):
        renderer = self.__set_indigo_renderer()
        common_molecules = []
        for entry in self.layout.meta_data['stock_solutions']:
            if 'all' in entry['wells']:
                common_molecules.append(entry['compound'])
        structures_buffers_dict = self.__create_structures_buffers_dict(renderer)
        width_factor, height_factor = self.__get_scaling(structures_buffers_dict)
        rows_table_list, columns_table_list = self.__load_images_for_tables(structures_buffers_dict, width_factor, height_factor)
        if common_molecules:
            common_mol_image = self.__create_common_molecules_image(common_molecules, width_factor, height_factor)
        else:
            common_mol_image = Paragraph('')
        rows_table = Table(rows_table_list)
        columns_table = Table(columns_table_list)
        return rows_table, columns_table, common_mol_image

    def __set_indigo_renderer(self):
        indigo.setOption("render-output-format", 'png')
        indigo.setOption('render-relative-thickness', 1.5)
        indigo.setOption('render-bond-length', 100)
        indigo.setOption('render-margins', 10, 10)
        return renderer

    def __create_structures_buffers_dict(self, renderer:IndigoRenderer):
        molecules = {'rows': self.layout.rows, 'columns': self.layout.columns}

        for key, item in molecules['rows'].items():
            mols = [item]
            for entry in self.layout.meta_data['stock_solutions']:
                if str(key) in entry['wells']:
                    if entry['compound'] not in mols:
                        mols.append(entry['compound'])
            molecules['rows'][key] = mols

        for key, item in molecules['columns'].items():
            mols = [item]
            for entry in self.layout.meta_data['stock_solutions']:
                if str(key) in entry['wells']:
                    if entry['compound'] not in mols:
                        mols.append(entry['compound'])
            molecules['columns'][key] = mols

        structures_buffers_dict = {'rows': {}, 'columns': {}}
        for key, item in molecules.items():
            for pos, smiles_list in item.items():
                if not isinstance(smiles_list, list):
                    smiles_list = [smiles_list]
                count = 0
                mol_array = indigo.createArray()
                strings = []
                if 'x' in smiles_list or 'X' in smiles_list:
                    pass
                else:
                    for smiles in smiles_list:
                        if '*' in smiles:
                            strings.append(smiles)
                        else:
                            mol = indigo.loadMolecule(smiles)
                            mol_array.arrayAdd(mol)
                            count += 1
                    buffer = renderer.renderGridToBuffer(mol_array, None, count)
                    if strings:
                        structures_buffers_dict[key][pos] = [io.BytesIO(buffer), strings]
                    else:
                        structures_buffers_dict[key][pos] = [io.BytesIO(buffer)]
        return structures_buffers_dict

    def __get_scaling(self, structures_buffers_dict:dict):
        width = self.doc.width/2-1*cm
        height = 2*cm
        buffers = list(structures_buffers_dict['rows'].values()) + list(structures_buffers_dict['columns'].values())
        buffers = [i[0] for i in buffers]
        width_factors = []
        height_factors = []
        for buffer in buffers:
            img = utils.ImageReader(buffer)
            iw, ih = img.getSize()
            width_factors.append(width/iw)
            height_factors.append(height/ih)
        width_factor = min(width_factors)
        height_factor = min(height_factors)
        return width_factor, height_factor

    def __scale_images(self, buffers:list[io.BytesIO], width_factor: float, height_factor: float):
        images = []
        for i, buffer in enumerate(buffers):
            img = utils.ImageReader(buffer)
            iw, ih = img.getSize()
            width_aspect = ih/float(iw)
            height_aspect = iw/float(ih)
            if width_factor <= height_factor:
                width = width_factor*iw
                height = (width*width_aspect)
            else:
                height = height_factor*ih
                width = (height*height_aspect)
            images.append(Image(buffer, width=width, height=height))
        return images

    def __load_images_for_tables(self, structures_buffers_dict:dict, width_factor:float, height_factor:float):
        lists = []
        for key in ['rows', 'columns']:
            buffers = structures_buffers_dict[key].values()
            buffers = [i[0] for i in buffers]
            images = self.__scale_images(buffers, width_factor, height_factor)
            for i, pos in enumerate(structures_buffers_dict[key]):
                structures_buffers_dict[key][pos][0] = images[i]
            lists.append(list(structures_buffers_dict[key].items()))
        return lists

    def __create_product_image(self, product:str):
        renderer = self.__set_indigo_renderer()
        mol = indigo.loadMolecule(product)
        buffer = io.BytesIO(renderer.renderToBuffer(mol))
        w, h = self.metadata.wrap(self.doc.width*2/3, self.doc.height)
        h = h - 25
        image = self.__scale_image_to_target(buffer, w, h)
        return image

    def __scale_image_to_target(self, buffer:io.BytesIO, target_width:float, target_height:float):
        img = utils.ImageReader(buffer)
        iw, ih = img.getSize()
        width_factor = target_width/iw
        height_factor = target_height/ih
        if width_factor <= height_factor:
            width_aspect = ih / float(iw)
            width = width_factor * iw
            height = (width * width_aspect)
        else:
            height_aspect = iw / float(ih)
            height = height_factor * ih
            width = (height * height_aspect)
        return Image(buffer, width=width, height=height)

    def __create_common_molecules_image(self, common_molecules: list, width_factor:float, height_factor:float):
        mol_array = indigo.createArray()
        for smiles in common_molecules:
            try:
                mol = indigo.loadMolecule(smiles)
                mol_array.arrayAdd(mol)
            except:
                print(f'{smiles} could not be rendered bei Indigo.')
        buffer = renderer.renderGridToBuffer(mol_array, None, len(common_molecules))
        buffer = io.BytesIO(buffer)
        image = self.__scale_image_by_factor(buffer, width_factor, height_factor)
        return image

    def __scale_image_by_factor(self, buffer:io.BytesIO, width_factor:float, height_factor:float):
        img = utils.ImageReader(buffer)
        iw, ih = img.getSize()
        if width_factor <= height_factor:
            width_aspect = ih / float(iw)
            width = width_factor * iw
            height = (width * width_aspect)
        else:
            height_aspect = iw / float(ih)
            height = height_factor * ih
            width = (height * height_aspect)
        return Image(buffer, width=width, height=height)

    def create_conditions_table(self):
        rows_list = []
        for key, value in self.layout.meta_data['conditions'].items():
            if key == 'vessel':
                row = [f'{key}:', f'{value["model"]}']
            else:
                row = [f'{key}:', f'{value}']
            rows_list.append(row)
        table = Table(rows_list, rowHeights=12, hAlign='LEFT')
        table.setStyle(TableStyle([('ALIGNMENT', (0, 0), (-1, -1), 'LEFT'),
                                ('FONTNAME', (0, 0), (-1, -1), 'Helvetica')]))
        return table

    def __create_analysis_table(self):
        rows_list = []
        for key in self.layout.meta_data['analysis'][0].keys():
            row = []
            for analysis in self.layout.meta_data['analysis']:
                row.extend([f'{key}:', f'{analysis[key]}'])
            rows_list.append(row)
        table = Table(rows_list, rowHeights=12, hAlign='LEFT')
        table.setStyle(TableStyle([('ALIGNMENT', (0, 0), (-1, -1), 'LEFT'),
                                   ('FONTNAME', (0, 0), (-1, -1), 'Helvetica')]))
        return table
