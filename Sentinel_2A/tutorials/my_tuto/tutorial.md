---
layout: tutorial_hands_on

title: Sentinel_2A_biodiversity
zenodo_link: ''
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- yvanlebras

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

This tutorial will guide you on getting Sentinel 2A data and precessing them in order to calculate and visualize biodiversity indicators.
We’ll be using data extracted from the Copernicus portal, Scihub.  Fisrt those data will be reformated. We’ll explore this dataset in the view of making biodiversity analyses so we will compute alpha and beta indicators, calculate Shannon, Hill and others indices and finally plot how the data is distributed through space and time, etc … .


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a name (example: “Sentinel 2A data for biodiversity tutorial”) for you to find it again later if needed.
>    {% snippet faqs/galaxy/histories_create_new.md box_type="none" %}
> 2. Download the files from [Scihub]({{ https://scihub.copernicus.eu/dhus/#/home }}) or from
>    [PEPS]({{ https://peps.cnes.fr/rocket/#/search?maxRecords=50&page=1 }} :
>
>    You will need to to create an account for either of these platform.
>    
>    You will download a zip folder. Keep it that way. 
> 3. Upload the zip folder
>
> ### {% icon tip %} Tip: Importing data from your computer
>
>    * Open the Galaxy Upload Manager {% icon galaxy-upload %}
>    * Select **Choose local files**
>    * Browse in your computer and get the downloaded zip folder
>    * Press **Start**
{: .tip}
>
> 4. Rename the datasets “sentinel_2a_data.zip" for example and preview your dataset
>
>
{: .hands_on}

# Reformating data

>   Unzips, reads and preprocesses Sentinel 2A zip folder data.

## Sub-step with **Preprocessing sentinel 2 data**

> ### {% icon hands_on %} Hands-on: Preprocess
>
> 1. {% tool [Preprocessing sentinel 2 data](SRS_preprocess_S2) %} with the following parameters:
>    - {% icon param-file %} *"Input data"*: `sentinel_2a_data.zip` (Input dataset)
>    - {% icon param-select %} *"Where does your data come from ?"*: 'From Scihub or Peps'
>
> 2. Click on **Execute** 
>
>    > ### {% icon comment %} Comment
>    >
>    > The interesting output is the ENVI image format which is a flat-binary raster file with an accompanying ASCII header file. The data are stored as a binary stream of bytes in a file without extension
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> >
> {: .solution}
>
{: .question}

# Producing biodiversity indicateurs

>   You can choose to compute spectral and biodiversity indicators either for global remote sensing data or for a canopy.

## Sub-step with **Compute spectral indices**

> ### {% icon hands_on %} Hands-on: spectral indices
>
> 1. {% tool [Compute spectral indices](SRS_canopy_indices) %} with the following parameters:
>    - {% icon param-file %} *"Input raster"*: `output_refl` (output of **Preprocessing sentinel 2 data** {% icon tool %})
>    - {% icon param-file %} *"Input raster header"*: `output_refl` (output of **Preprocessing sentinel 2 data** {% icon tool %})
>    - {% icon param-select %} *"Input the type of indice you want"*: 'NDVI'
>
>    ***TODO***: *Check that the "Input raster" is bil datatype and that "Input raster header" is a hdr datatype*
> ### {% icon tip %} Tip: Checking datatype
>
>    * Go on your raster data 
>    * Click on {% icon galaxy-pencil %} to edit it
>    * Click on {% icon galaxy-chart-select-data %} Datatypes
>    * On "New Type" **Select** bil 
>    * Press **Save**
>
>    > ### {% icon comment %} Comment
>    >
>    > Do the same for the raster header with the datatype hdr
>    {: .comment}
{: .tip}
>
>    > ### {% icon comment %} Comment
>    >
>    > You can choose whichever indice you want
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> >
> {: .solution}
>
{: .question}

## Sub-step with **Compute biodiversity indices**

> ### {% icon hands_on %} Hands-on: Biodiversity indicators for global remote sensing data
>
> 1. {% tool [Compute biodiversity indices](SRS_global_indices) %} with the following parameters:
>    - {% icon param-file %} *"Input raster"*: `output_refl` (output of **Preprocessing sentinel 2 data** {% icon tool %})
>    - {% icon param-file %} *"Input raster header"*: `output_refl` (output of **Preprocessing sentinel 2 data** {% icon tool %})
>    - {% icon param-text %} *"Write a number of the value of alpha"*: '1'
>
>
>    > ### {% icon comment %} Comment
>    >
>    > Same as the compute spectral indices make sure you have the right datatypes bil and hdr.
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> >
> {: .solution}
>
{: .question}

## Sub-step with **Processing remote sensing data**

> ### {% icon hands_on %} Hands-on: Biodiversity indicators for canopy remote sensing data
>
> 1. {% tool [Processing remote sensing data](SRS_process_data) %} with the following parameters:
>    - {% icon param-file %} *"Input raster"*: `output_refl` (output of **Preprocessing sentinel 2 data** {% icon tool %})
>    - {% icon param-file %} *"Input raster header"*: `output_refl` (output of **Preprocessing sentinel 2 data** {% icon tool %})
>    - {% icon param-file %} *"Plots folder zip"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> >
> {: .solution}
>
{: .question}


# Conclusion
{:.no_toc}

Sentinel 2A biodiversity data tutorial finished !
