const fs = require('fs');
const path = require('path');
const puppeteer = require('puppeteer');

// Define the directories for HTML input and PDF output
const inputDir = '/Users/casper/Documents/data/webapp-paper/html_one_pagers/';
const outputDir = '/Users/casper/Documents/data/webapp-paper/pdf_one_pagers/';

// Function to generate PDF from HTML
async function generatePDF(htmlFilePath, pdfFilePath) {
    const browser = await puppeteer.launch();
    const page = await browser.newPage();
    
    // Load HTML content from a file
    const html = fs.readFileSync(htmlFilePath, 'utf8');
    
    // Set the HTML content on the page
    await page.setContent(html, { waitUntil: 'networkidle0' });

    // Generate the PDF and save it
    await page.pdf({
        path: pdfFilePath,  // Output file path
        format: 'A4',  // Page format
        printBackground: true  // Include background images and colors
    });

    await browser.close();
}

// Read all HTML files from the input directory
fs.readdir(inputDir, async (err, files) => {
    if (err) {
        console.error('Error reading input directory:', err);
        return;
    }

    // Filter HTML files
    const htmlFiles = files.filter(file => file.endsWith('.html'));

    // Loop through each HTML file
    for (const file of htmlFiles) {
        const htmlFilePath = path.join(inputDir, file);
        const pdfFilePath = path.join(outputDir, file.replace('.html', '.pdf'));

        try {
            console.log(`Generating PDF for: ${file}`);
            await generatePDF(htmlFilePath, pdfFilePath);
            console.log(`PDF saved: ${pdfFilePath}`);
        } catch (error) {
            console.error(`Error generating PDF for ${file}:`, error);
        }
    }
});