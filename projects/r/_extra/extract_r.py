import re

with open('/Users/jay/Desktop/websiteV02/projects/r/airbnb_sf/AIRBNB_SanFran.qmd', 'r') as f:
    content = f.read()

# Find all R code chunks: match ```{r ...} ... ```
chunks = re.findall(r'```\{r.*?\}\n(.*?)```', content, re.DOTALL)

with open('/Users/jay/Desktop/websiteV02/projects/r/airbnb_sf/temp_airbnb.R', 'w') as f:
    for i, chunk in enumerate(chunks):
        f.write(f"# Chunk {i+1}\n")
        f.write(chunk)
        f.write("\n")
