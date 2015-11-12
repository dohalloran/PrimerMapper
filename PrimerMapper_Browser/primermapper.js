window.onload = function () {
    var fileInput = document.getElementById('fileInput');
    var fileDisplayArea = document.getElementById('fileDisplayArea');

    fileInput.addEventListener('change', function (e) {
        var file = fileInput.files[0];

        if (file) {
            var reader = new FileReader();

            reader.onload = function (e) {
                fileDisplayArea.innerText = reader.result;
            }

            reader.readAsText(file);
        } else {
            fileDisplayArea.innerText = "File not supported!";
        }
    });
}


var nn_s = {
    'AA': '240',
    'AC': '173',
    'AG': '208',
    'AT': '239',
    'AN': '215',
    'CA': '129',
    'CC': '266',
    'CG': '278',
    'CT': '208',
    'CN': '220',
    'GA': '135',
    'GC': '267',
    'GG': '266',
    'GT': '173',
    'GN': '210',
    'TA': '169',
    'TC': '135',
    'TG': '129',
    'TT': '240',
    'TN': '168',
    'NA': '168',
    'NC': '210',
    'NG': '220',
    'NT': '215',
    'NN': '203',
    'aa': '240',
    'ac': '173',
    'ag': '208',
    'at': '239',
    'an': '215',
    'ca': '129',
    'cc': '266',
    'cg': '278',
    'ct': '208',
    'cn': '220',
    'ga': '135',
    'gc': '267',
    'gg': '266',
    'gt': '173',
    'gn': '210',
    'ta': '169',
    'tc': '135',
    'tg': '129',
    'tt': '240',
    'tn': '168',
    'na': '168',
    'nc': '210',
    'ng': '220',
    'nt': '215',
    'nn': '203'
}
var nn_h = {
    'AA': '91',
    'AC': '65',
    'AG': '78',
    'AT': '86',
    'AN': '80',
    'CA': '58',
    'CC': '110',
    'CG': '119',
    'CT': '78',
    'CN': '91',
    'GA': '56',
    'GC': '111',
    'GG': '110',
    'GT': '65',
    'GN': '85',
    'TA': '60',
    'TC': '56',
    'TG': '58',
    'TT': '91',
    'TN': '66',
    'NA': '66',
    'NC': '85',
    'NG': '91',
    'NT': '80',
    'NN': '80',
    'aa': '91',
    'ac': '65',
    'ag': '78',
    'at': '86',
    'an': '80',
    'ca': '58',
    'cc': '110',
    'cg': '119',
    'ct': '78',
    'cn': '91',
    'ga': '56',
    'gc': '111',
    'gg': '110',
    'gt': '65',
    'gn': '85',
    'ta': '60',
    'tc': '56',
    'tg': '58',
    'tt': '91',
    'tn': '66',
    'na': '66',
    'nc': '85',
    'ng': '91',
    'nt': '80',
    'nn': '80'
}


function loaderMy() {

    var gene_length = "";
    var start_F = "";
    var start_R = "";
    var length_F = "";
    var length_R = "";
    var dnaForPrimer_with_label = "";
    var dnaSequence = [];
    var splittedLines = document.getElementById("page-wrapper").innerText;
    dnaSequence = splittedLines.split('\n');

    for (var i = 0; i < dnaSequence.length; i++) {
        if (dnaSequence[i].match(/#/)) {
            gene_length = dnaSequence[i];
        } else if (dnaSequence[i].match(/FS/)) {
            start_F = dnaSequence[i];
        } else if (dnaSequence[i].match(/RS/)) {
            start_R = dnaSequence[i];
        } else if (dnaSequence[i].match(/LF/)) {
            length_F = dnaSequence[i];
        } else if (dnaSequence[i].match(/LR/)) {
            length_R = dnaSequence[i];
        } else if (dnaSequence[i].match(/>/)) {
            gene_name = dnaSequence[i];
        } else if (dnaSequence[i].match(/SEQ/)) {
            dnaForPrimer_with_label = dnaSequence[i];
        }
    }


    var gene_len = gene_length.replace("#", "");
    var geneName = gene_name.replace(">", "");

    var dnaForPrimers = dnaForPrimer_with_label.replace("SEQ", "");
    dnaForPrimers = dnaForPrimers.replace(/\s+/g, '');

    var replaced_start_F = start_F.replace("FS", "");
    var replaced_start_R = start_R.replace("RS", "");
    var replaced_length_F = length_F.replace("LF", "");
    var replaced_length_R = length_R.replace("LR", "");

    var primer_start_F = replaced_start_F.split(',');
    var primer_length_F = replaced_length_F.split(',');
    var primer_start_R = replaced_start_R.split(',');
    var primer_length_R = replaced_length_R.split(',');



    //////////////////////////////

    /////DRAW RULER 


    var canvas = new fabric.Canvas('canvas2');

    line_length = gene_len;

    adjusted_length = (line_length / 666) * 167;

    canvas.add(new fabric.Line([0, 0, adjusted_length, 0], {
        left: 50,
        top: 0,
        stroke: '#d89300',
        strokeWidth: 3
    }));

    $('#canvas_container').css('overflow-x', 'scroll');
    $('#canvas_container').css('overflow-y', 'hidden');

    drawRuler();


    function drawRuler() {
        $("#ruler[data-items]").val(line_length / 200);

        $("#ruler[data-items]").each(function () {
            var ruler = $(this).empty(),
                len = Number($("#ruler[data-items]").val()) || 0,
                item = $(document.createElement("li")),
                i;
            ruler.append(item.clone().text(1));
            for (i = 1; i < len; i++) {
                ruler.append(item.clone().text(i * 200));
            }
            ruler.append(item.clone().text(i * 200));
        });
    }


    //////////////////////////////

    /////MAP MOUSE POSITION

    function mousePos() {
        canvas.add(new fabric.Text('base pairs', {
            fontStyle: 'italic',
            fontFamily: 'Hoefler Text',
            fontSize: 12,
            left: 0,
            top: 10,
            hasControls: false,
            selectable: false
        }));

        canvas.on('mouse:move', function (options) {
            getMouse(options);
        });


        function getMouse(options) {
            canvas.getObjects('text')[0].text =
                "base pairs: " + Math.round(Math.round(funcAdjust(options.e.clientX) - 13.5) * 15.9);
            canvas.renderAll();
        }

    }

    mousePos();

    //////////////////////////////
    /////DRAW PRIMERS 

    var canvasNewer = new fabric.Canvas('canvas3', {
        hoverCursor: 'pointer',
        selection: false
    });

    fabric.Object.prototype.transparentCorners = false;

    function drawMyPrimers() {
        canvasNewer.add(new fabric.Text('draw primers:', {
            fontStyle: 'italic',
            fontFamily: 'Hoefler Text',
            fontSize: 12,
            left: 0,
            top: 0,
            hasControls: false,
            selectable: false
        }));

        canvasNewer.getObjects('text')[0].text = "draw primers:";

        var getPosit = [];
        canvasNewer.on('mouse:down', function (o) {
            isDown = true;
            var pointer = canvasNewer.getPointer(o.e);
            var posX1 = pointer.x;
            getPosit.push(posX1);
            var points = [pointer.x, pointer.y, pointer.x, pointer.y];
            line = new fabric.Line(points, {
                strokeWidth: 2,
                fill: 'green',
                stroke: 'green',
                originX: 'center',
                originY: 'center'
            });
            canvasNewer.add(line);
        });


        canvasNewer.on('mouse:move', function (o) {
            if (!isDown) return;
            var pointer = canvasNewer.getPointer(o.e);
            line.set({
                x2: pointer.x,
                y2: pointer.y
            });
            canvasNewer.renderAll();
        });

        canvasNewer.on('mouse:up', function (o) {
            var pointer = canvasNewer.getPointer(o.e);
            var posX2 = pointer.x;
            getPosit.push(posX2);
            isDown = false;
            var startPrimer = getPosit[0];
            var endPrimer = getPosit[1];
            var startScaled = Math.round(Math.round(funcAdjust(startPrimer) - 12.5) * 15.9);
            var endScaled = Math.round(Math.round(funcAdjust(endPrimer) - 12.5) * 15.9);
            var primer = dnaForPrimers.substring(startScaled, endScaled);

            if (primer.length > 36) {
            if (endScaled >= startScaled) {
                alert("Forward Primer: " + primer + "\n" + "GC percent: " + gcPercent(primer) + "%" + "\n" + "Primer Tm: " + calcTmLong(primer) + "\u2103");
            } else {
                alert("Reverse Primer: " + reverse(primer) + "\n" + "GC percent: " + gcPercent(reverse(primer)) + "%" + "\n" + "Primer Tm: " + calcTmLong(reverse(primer)) + "\u2103");
            }
            } else {
               if (endScaled >= startScaled) {
                alert("Forward Primer: " + primer + "\n" + "GC percent: " + gcPercent(primer) + "%" + "\n" + "Primer Tm: " + calcTm(primer) + "\u2103");
            } else {
                alert("Reverse Primer: " + reverse(primer) + "\n" + "GC percent: " + gcPercent(reverse(primer)) + "%" + "\n" + "Primer Tm: " + calcTm(reverse(primer)) + "\u2103");
            }     
            }
            getPosit = [];
        });
    }

    drawMyPrimers();

    //////////////////////////////


    /////MAP PRIMERS

    var canvasNew = document.getElementById('canvas');
    var context = canvasNew.getContext('2d');


    context.font = "12px Georgia";
    context.fillText("PrimerMapper results for " + geneName, 10, 20);

    function Line(x1, y1, x2, y2) {
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
    }
    Line.prototype.drawWithArrowheads = function (ctx) {

        ctx.strokeStyle = "blue";
        ctx.fillStyle = "blue";
        ctx.lineWidth = 1;

        ctx.beginPath();
        ctx.moveTo(this.x1, this.y1);
        ctx.lineTo(this.x2, this.y2);
        ctx.stroke();

        var endRadians = Math.atan((this.y2 - this.y1) / (this.x2 - this.x1));
        endRadians += ((this.x2 > this.x1) ? 90 : -90) * Math.PI / 180;
        this.drawArrowhead(ctx, this.x2, this.y2, endRadians);

    }
    Line.prototype.drawArrowhead = function (ctx, x, y, radians) {
        ctx.save();
        ctx.beginPath();
        ctx.translate(x, y);
        ctx.rotate(radians);
        ctx.moveTo(0, 0);
        ctx.lineTo(3, 10);
        ctx.lineTo(-3, 10);
        ctx.closePath();
        ctx.restore();
        ctx.fill();
    }


    function LineR(x1, y1, x2, y2) {
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
    }

    LineR.prototype.drawWithArrowheads = function (ctz) {

        ctz.strokeStyle = "red";
        ctz.fillStyle = "red";
        ctz.lineWidth = 1;

        ctz.beginPath();
        ctz.moveTo(this.x1, this.y1);
        ctz.lineTo(this.x2, this.y2);
        ctz.stroke();

        var startRadians = Math.atan((this.y2 - this.y1) / (this.x2 - this.x1));
        startRadians += ((this.x2 > this.x1) ? -90 : 90) * Math.PI / 180;
        this.drawArrowhead(ctz, this.x1, this.y1, startRadians);

    }
    LineR.prototype.drawArrowhead = function (ctz, x, y, radians) {
        ctz.save();
        ctz.beginPath();
        ctz.translate(x, y);
        ctz.rotate(radians);
        ctz.moveTo(0, 0);
        ctz.lineTo(3, 10);
        ctz.lineTo(-3, 10);
        ctz.closePath();
        ctz.restore();
        ctz.fill();
    }


    if (document.getElementById('rad1').checked) {

        var k;
        var leftOver_F = 25;
        for (k = 0; k < primer_start_F.length; k++) {
            var adjusted_primer_start_F = funcAdjust(primer_start_F[k]);
            var adjusted_primer_length_F = funcAdjust(primer_length_F[k]);
            leftOver_F += 15;

            var line = new Line(30 + adjusted_primer_start_F, leftOver_F, 40 + adjusted_primer_start_F + adjusted_primer_length_F, leftOver_F);
            line.drawWithArrowheads(context);
            context.font = "8px Arial";
            context.fillText(primer_start_F[k], adjusted_primer_start_F + 5, leftOver_F + 3);

        }

        var j;
        var leftOver_R = 25;
        for (j = 0; j < primer_start_R.length; j++) {
            var adjusted_primer_start_R = funcAdjust(primer_start_R[j]);
            var adjusted_primer_length_R = funcAdjust(primer_length_R[j]);
            leftOver_R += 15;

            var line2 = new LineR(30 + adjusted_primer_start_R, leftOver_R, 40 + adjusted_primer_start_R + adjusted_primer_length_R, leftOver_R);
            line2.drawWithArrowheads(context);
            context.font = "8px Arial";
            context.fillText(primer_start_R[j], adjusted_primer_start_R + 62, leftOver_R + 3);


        }
    } else if (document.getElementById('rad2').checked) {

        var k;
        var leftOver_F = 25;
        for (k = 0; k < primer_start_F.length; k++) {
            var adjusted_primer_start_F = funcAdjust(primer_start_F[k]);
            var adjusted_primer_length_F = funcAdjust(primer_length_F[k]);
            leftOver_F += 15;

            context.beginPath();
            context.moveTo(30 + adjusted_primer_start_F, leftOver_F);
            context.lineTo(30 + adjusted_primer_start_F + adjusted_primer_length_F, leftOver_F);
            context.lineWidth = 5;
            context.strokeStyle = '#0000FF';
            context.stroke();
            context.font = "8px Arial";
            context.fillText(primer_start_F[k], adjusted_primer_start_F + 5, leftOver_F + 3);


        }

        var j;
        var leftOver_R = 25;
        for (j = 0; j < primer_start_R.length; j++) {
            var adjusted_primer_start_R = funcAdjust(primer_start_R[j]);
            var adjusted_primer_length_R = funcAdjust(primer_length_R[j]);
            leftOver_R += 15;

            context.beginPath();
            context.moveTo(30 + adjusted_primer_start_R, leftOver_R);
            context.lineTo(30 + adjusted_primer_start_R + adjusted_primer_length_R, leftOver_R);
            context.lineWidth = 5;
            context.strokeStyle = '#FF0000';
            context.stroke();
            context.font = "8px Arial";
            context.fillText(primer_start_R[j], adjusted_primer_start_R + 62, leftOver_R + 3);

        }

    }

    //Calculate Tm of Primer if >36bps
    function calcTmLong(dna) {
        var numberTmLong = 81.5 + (16.6 * (Math.log(50 / 1000.0) / Math.log(10))) +
            (41.0 * (gcPercent(dna) / 100)) - (600.0 / dna.length);
        return numberTmLong.toFixed(2);
    }

    //Calculate Tm of Primer if <=36bps
    function calcTm(sequence) {
        var dH = 0;
        var dS = 108;
        var i;
        // Compute dH and dS
        for (i = 0; i < (sequence.length - 1); i++) {
            var pair = sequence.substr(i, 2);
            dH += parseInt(nn_h[pair], 10);
            dS += parseInt(nn_s[pair], 10);
        }
        dH *= -100.0;
        dS *= -0.1;
        var numberTm = dH / (dS + 1.987 * Math.log(100 / 4000000000.0)) - 273.15 +
            16.6 * (Math.log(50 / 1000.0) / Math.log(10));
        return numberTm.toFixed(2);
    }


    //Calculate GC Percentage of Primer
    function gcPercent(dna) {
        var Arr_Primer = dna.split("");
        var gCount = 0;
        var tCount = 0;
        var cCount = 0;
        var aCount = 0;
        for (var i = 0; i < Arr_Primer.length; i++) {
            if (Arr_Primer[i] === 'A') {
                aCount++;
            } else if (Arr_Primer[i] === 'C') {
                cCount++;
            } else if (Arr_Primer[i] === 'T') {
                tCount++;
            } else if (Arr_Primer[i] === 'G') {
                gCount++;
            }
        }
        var number = ((gCount + cCount) / Arr_Primer.length) * 100;
        return number.toFixed(2);
    }


    //Reverse Complement Primer
    function reverse(s) {
        return s.split("").map(complement).reverse().join("");
    }


    function complement(nucleotide) {
        var complements = {
            'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A'
        };

        return complements[nucleotide];
    }


    function funcAdjust(start_bp) {
        var adjusted = (start_bp / 666) * 167;
        return adjusted;
    }

    ///////////
    //////////magnifying glass
    if (document.getElementById('lens').checked) {
        var zoom = document.getElementById("zoom");
        var zoomCtx = zoom.getContext("2d");
        canvasNew.addEventListener("mousemove", function (e) {
            console.log(e);
            //zoomCtx.fillStyle = "white";
            zoomCtx.clearRect(0, 0, zoom.width, zoom.height);
            zoomCtx.fillStyle = "transparent";
            //zoomCtx.fillRect(0,0, zoom.width, zoom.height);
            zoomCtx.drawImage(canvasNew, e.x, e.y - 120, 1600, 800, 0, 0, 3200, 1600);
            console.log(zoom.style);
            zoom.style.top = e.pageY - 5 + "px"
            zoom.style.left = e.pageX + 20 + "px"
            zoom.style.display = "block";
        });

        canvasNew.addEventListener("mouseout", function () {
            zoom.style.display = "none";
        });
    }
    //////////
    //////////
    ///////////


}