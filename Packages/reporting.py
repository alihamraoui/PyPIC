from datetime import datetime
import pandas as pd

def reporting (table) :
    """
    This function product html repport with results of interaction
    """
    #Converting to HTML# saving plot image to local file
    #image = plot.get_figure().savefig('plot.png')

    #image_tag = '<img src="plot.png">'#writing HTML Content

    infile = '../test/results/' + table
    df = pd.read_table (infile)

    heading = '<h1> PyPIC Protein Interaction Calculator Report </h1>'

    subheading = '<h3> Hydrophobic interactions </h3>'# Using .now() from datetime library to add Time stamp

    paragraph = '<p> The following residues are considered to participate in interactions if they fall within 5Å range.  \n' + \
                'ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, TYR.  \n' + \
                'Reference: Kyte and Doolittle. </p>'

    now = datetime.now()

    current_time = now.strftime("%m/%d/%Y %H:%M:%S")

    header = '<div class="top">' + heading + subheading +'</div>'

    footer = '<div class="bottom"> <h3> This Report has been Generated on '+ current_time +'</h3> </div>'

    content = paragraph + '<div class="table"> '+df.to_html( index=False)+' </div> \n <div class="chart"> '+ 'image_tag' +'</div>'# Concating everything to a single string

    html = header + content + footer

    outfile = '../test/results/report.html'
    with open(outfile,'w+') as file:
        file.write(html)