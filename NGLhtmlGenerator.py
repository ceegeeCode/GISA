# -*- coding: utf-8 -*-
"""
Created on Sat Feb 01 20:28:52 2020

@author: Christian Grønbæk
"""


'''
Script for generating HTML-code for rendering interactive 3d-plots, all based
on the NGL viewer (A.Rose, A.Bradley, Y.Valasatava, J.Duarte A.Prlić, P.Rose, "NGL Viewer: Web-based molecular graphics for large complexes", Bioinformatics, 2018
and A.Rose, P.Hildebrand, "NGL Viewer: a web application for molecular visualization", Nucleic Acids Research, 2015).



Usage:

import NGLhtmlGenerator as nglGenerator

rootOutput = r'C:\Users\Christian\Sync\Bioinformatics\papers\geometric_anomalies_in_folds\plots\htmlsForNGL'

pdbId = '1dif'
chainId = ''
segment1 = [10,20]
segment2 = [30,50]
fileName = r'\\' + pdbId + '.htm'
fileNameOutput = rootOutput +  fileName

html = nglGenerator.generateHTML(pdbId = pdbId, chainId = chainId, segment1 = segment1, segment2 = segment2, rootOutput = rootOutput)

The html-code generated is now stored in the output file (rootOutput + pdbId + chainId + 'htm') ; opening it in an internet browser will produce a cartoon 3d-plot of the 
chain/structure having the two segments high-lighted in blue/red (the N-terminus closes being blue).

'''

def generateHTML(pdbId, chainId, segment1, segment2, strValue, rootOutput):
    '''Function returning html-code by which a NGL viewer 3d-plot 
    can be generated; the html can be opened in a internet-browser.
    
    strValue: typically the mutual writhe of the subchain pair (given by the segments).
    '''
    
    if chainId == '>':
        chainId = ''
    
    seg1String = str(segment1[0]) + '-' + str(segment1[1]) 
    seg2String = str(segment2[0]) + '-' + str(segment2[1]) 
    
    htmlString = '<!DOCTYPE html>' + '\n'
    htmlString += '<html lang="en">' + '\n'
    htmlString += '<head>' + '\n'
    htmlString += '  <meta charset="utf-8">' + '\n'
    htmlString += '</head>' + '\n'
    htmlString += '<body>' + '\n'
    htmlString += '  <pre>' + '\n'
    htmlString += '               ' + pdbId + chainId + ';' + seg1String + ';' + seg2String + '; w:' + strValue  +'\n'
    htmlString += '</pre>' + '\n'
    htmlString += '  <script src="https://unpkg.com/ngl"></script>' + '\n'
    htmlString += '  <script>' + '\n'
    htmlString += '    document.addEventListener("DOMContentLoaded", function () {' + '\n'
    htmlString += '            var stage = new NGL.Stage("viewport",  {backgroundColor:"white"});' + '\n'
    htmlString += '            var schemeId = NGL.ColormakerRegistry.addSelectionScheme([' + '\n'
    htmlString += '  ["blue", "' + seg1String + '"],' + '\n'
    htmlString += '  ["red", "' + seg2String + '"], ' + '\n'
    htmlString += '  ["khaki", "*"]' + '\n'
    htmlString += '  ]);' + '\n'
    htmlString += '  ' + '\n'
    htmlString += 'stage.loadFile("rcsb://' + pdbId + '").then(function (o) {' + '\n'
    htmlString += '  o.addRepresentation("cartoon",  {sele: ":'+ chainId +'", color: schemeId});  // pass schemeId here' + '\n'
    htmlString += '  o.autoView();' + '\n'
    htmlString += '});' + '\n'
    htmlString += '' + '\n'
    htmlString += '' + '\n'
    htmlString += '});' + '\n'
    htmlString += '' + '\n'
    htmlString += '' + '\n'
    htmlString += '  </script>' + '\n'
    htmlString += '  <div id="viewport" style="width:600px; height:500px;"></div>' + '\n'
    htmlString += '</body>' + '\n'
    htmlString += '</html>' + '\n'
    
    fileName = pdbId + chainId + '_' + seg1String + '_' + seg2String + '.htm'
    fileNameOutput = rootOutput + r'\\' + fileName
    print fileNameOutput
    outFile = open(fileNameOutput, 'w')
    outFile.write(htmlString)
    outFile.close()
    
    return htmlString
    
    