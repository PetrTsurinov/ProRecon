# ProRecon Java library
## Overview
ProRecon is a lightweight Java library that is designed to RECONstruct PROteins consider different mutations. <br />
Detailed presentation can be found here: <br />
Pdf - https://drive.google.com/open?id=0B1g2jSBV2uqrZWljM3c4S1oxbGc <br />
Video - https://drive.google.com/open?id=0B1V1JjzgsVlLTjk4UHlrcHl3NmM <br />

## How to Use
For now example can be found in temporary Main class or you can use jar with dependencies: <br />
java -cp prorecon.jar com.epam.prorecon.Main "dmel-all-chromosome-r606.fasta" "agnX1.model.2.snp-indels.vcf" "dmel-all-r6.06.LIMK1.gtf" "X" 12584385 12592193 <br />
where "dmel-all-chromosome-r606.fasta" is file with original sequence, <br />
"agnX1.model.2.snp-indels.vcf" is file with mutations list, <br />
"dmel-all-r6.06.LIMK1.gtf" is file with exons positions indexes, <br />
"X" is chromosome name, <br />
12584385 is DNA sequence start index, <br />
12592193 is DNA sequence end index.

## License
ProRecon is free software: You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 3 of the License.

This program is distributed WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should receive a copy of the GNU General Public License along with this program. If you did not, please see http://www.gnu.org/licenses/.
